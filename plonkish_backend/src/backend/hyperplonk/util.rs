use crate::{
    backend::{
        hyperplonk::{
            preprocessor::{compose, permutation_polys},
            prover::{
                instance_polys, lookup_compressed_polys, lookup_h_polys, lookup_m_polys,
                permutation_z_polys,
            },
            // anemoi::anemoi_hash_circuit_info,
        },
        mock::MockCircuit,
        PlonkishCircuit, PlonkishCircuitInfo,
    },
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{powers, PrimeField},
        chain,
        expression::{rotate::Rotatable, Expression, Query, Rotation},
        test::{rand_array, rand_idx, rand_vec},
        Itertools,
    },
    anemoi_hash::{AnemoiJive256, AnemoiJive, ApplicableMDSMatrix, MDSMatrix},
};


use halo2_curves::bn256::Fr;
use halo2_curves::ff::Field;
use num_integer::Integer;
use rand::RngCore;
use std::{
    array,
    collections::{HashMap, HashSet},
    hash::Hash,
    iter, mem,
};



/// Create a PlonkishCircuitInfo for Anemoi hash with TurboPlonk constraints
/// Based on the paper "An efficient verifiable state for zk-EVM and beyond from the Anemoi hash function"
pub fn anemoi_hash_circuit_info<F: PrimeField>(
    num_vars: usize,
    _num_rounds: usize,
) -> PlonkishCircuitInfo<F> {
    // We need 5 wires as described in the paper: w1, w2, w3, w4, wo
    // And 14 selector polynomials for the Anemoi constraints
    
    // Create expressions for wire polynomials
    let [w1, w2, w3, w4, wo] = &array::from_fn(|i| {
        Query::new(i, Rotation::cur())
    }).map(Expression::Polynomial);
    
    // Create expressions for selector polynomials
    let [q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4] = 
        &array::from_fn(|i| Query::new(i + 5, Rotation::cur())).map(Expression::Polynomial);
    
    // Create expressions for next rotation (for Anemoi constraints)
    let [w1_next, w2_next, w3_next, _wo_next] = &[
        Query::new(0, Rotation::next()),
        Query::new(1, Rotation::next()),
        Query::new(2, Rotation::next()),
        Query::new(4, Rotation::next()),
    ].map(Expression::Polynomial);
    
    // Base constraint (TurboPlonk gate)
    let base_constraint = q1 * w1 + q2 * w2 + q3 * w3 + q4 * w4
        + qm1 * w1 * w2 + qm2 * w3 * w4 + qc
        + qecc * w1 * w2 * w3 * w4 * wo
        - qo * wo;
    
    // Boolean constraints
    let bool_constraints = vec![
        qb * w2 * (w2 - Expression::one()),
        qb * w3 * (w3 - Expression::one()),
        qb * w4 * (w4 - Expression::one()),
    ];
    
    // Anemoi constraints (from Section 6 of the paper)
    let g = Expression::<F>::Constant(F::from(5u64)); // Generator
    let g_inv = Expression::<F>::Constant(F::from(5u64).invert().unwrap()); // Delta
    
    // Helper expressions for Anemoi round
    let c_prime_1 = w1 + w4 + g.clone() * (w2 + w3) + qprk3;
    let c_prime_2 = g.clone() * (w1 + w4) + (g.clone() * g.clone() + Expression::one()) * (w2 + w3) + qprk4;
    
    // First Anemoi equation
    let anemoi_1 = qprk3.clone() * (
        (c_prime_1.clone() - w3_next).pow(5)
        + g.clone() * c_prime_1.clone().pow(2)
        - (Expression::Constant(F::from(2u64)) * w1 + w4 
          + g.clone() * (Expression::Constant(F::from(2u64)) * w2 + w3) + qprk1)
    );
    
    // Second Anemoi equation  
    let anemoi_2 = qprk3.clone() * (
        (c_prime_2.clone() - wo).pow(5)
        + g.clone() * c_prime_2.clone().pow(2)
        - (g.clone() * (Expression::Constant(F::from(2u64)) * w1 + w4) 
           + (g.clone() * g.clone() + Expression::one()) * (Expression::Constant(F::from(2u64)) * w2 + w3) + qprk2)
    );
    
    // Third Anemoi equation
    let anemoi_3 = qprk3.clone() * (
        (c_prime_1 - w3_next.clone()).pow(5)
        + g.clone() * w3_next.clone().pow(2)
        + g_inv.clone()
        - w1_next
    );
    
    // Fourth Anemoi equation
    let anemoi_4 = qprk3.clone() * (
        (c_prime_2 - wo.clone()).pow(5)
        + g * wo.clone().pow(2)
        + g_inv
        - w2_next
    );
    
    // Collect all constraints
    let mut constraints = vec![base_constraint];
    constraints.extend(bool_constraints);
    constraints.extend(vec![anemoi_1, anemoi_2, anemoi_3, anemoi_4]);
    
    // Create preprocessed polynomials (selectors) - 14 selectors
    let preprocess_polys = vec![vec![F::ZERO; 1 << num_vars]; 14];
    
    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![0],
        preprocess_polys,
        num_witness_polys: vec![5], // w1, w2, w3, w4, wo
        num_challenges: vec![0],
        constraints,
        lookups: Vec::new(),
        permutations: Vec::new(),
        max_degree: Some(7), // Due to pow(5) in Anemoi constraints
    }
}

/// Generate a complete Anemoi hash circuit implementing the paper's specifications
/// Automatically determines the required num_vars based on gate count
pub fn rand_anemoi_hash_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    mut _preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    const NUM_ROUNDS: usize = 14; // From AnemoiJive256 specification
    const GATES_NEEDED: usize = NUM_ROUNDS + 3; // Input + rounds + output + padding
    
    // Calculate minimum num_vars needed
    let num_vars = (GATES_NEEDED as f64).log2().ceil() as usize;
    let size = 1 << num_vars;
    let usable_indices = R::from(num_vars).usable_indices();
    
    // Initialize polynomials
    let mut w1_values = vec![F::ZERO; size];
    let mut w2_values = vec![F::ZERO; size];
    let mut w3_values = vec![F::ZERO; size];
    let mut w4_values = vec![F::ZERO; size];
    let mut wo_values = vec![F::ZERO; size];
    
    // Selector polynomials (14 total)
    let mut q1 = vec![F::ZERO; size];
    let mut q2 = vec![F::ZERO; size];
    let mut q3 = vec![F::ZERO; size];
    let mut q4 = vec![F::ZERO; size];
    let mut qo = vec![F::ZERO; size];
    let mut qm1 = vec![F::ZERO; size];
    let mut qm2 = vec![F::ZERO; size];
    let mut qc = vec![F::ZERO; size];
    let mut qecc = vec![F::ZERO; size];
    let mut qb = vec![F::ZERO; size];
    let mut qprk1 = vec![F::ZERO; size];
    let mut qprk2 = vec![F::ZERO; size];
    let mut qprk3 = vec![F::ZERO; size];
    let mut qprk4 = vec![F::ZERO; size];
    
    let mut permutation = Permutation::default();
    
    // Use real Anemoi constants
    let generator = F::from(5u64);
    let generator_inv = F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap_or(generator.invert().unwrap());
    let generator_square_plus_one = F::from(26u64);
    
    // Real ROUND_KEYS_X and ROUND_KEYS_Y from AnemoiJive254
    let round_keys_x_strings = [
        ["37", "3751828524803055471428227881618625174556947755988347881191159153764975591158"],
        ["13352247125433170118601974521234241686699252132838635793584252509352796067497", "21001839722121566863419881512791069124083822968210421491151340238400176843969"],
        ["8959866518978803666083663798535154543742217570455117599799616562379347639707", "21722442537234642741320951134727484119993387379465291657407115605240150584902"],
        ["3222831896788299315979047232033900743869692917288857580060845801753443388885", "5574110054747610058729632355948568604793546392090976147435879266833412620404"],
        ["11437915391085696126542499325791687418764799800375359697173212755436799377493", "19347108854758320361854968987183753113398822331033233961719129079198795045322"],
        ["14725846076402186085242174266911981167870784841637418717042290211288365715997", "17733032409684964025894538244134113560864261458948810209753406163729963104066"],
        ["3625896738440557179745980526949999799504652863693655156640745358188128872126", "16641102106808059030810525726117803887885616319153331237086309361060282564245"],
        ["463291105983501380924034618222275689104775247665779333141206049632645736639", "9245970744804222215259369270991414441925747897718226734085751033703871913242"],
        ["17443852951621246980363565040958781632244400021738903729528591709655537559937", "18243401795478654990110719981452738859015913555820749188627866268359980949315"],
        ["10761214205488034344706216213805155745482379858424137060372633423069634639664", "18200337361605220875540054729693479452916227111908726624753615870884702413869"],
        ["1555059412520168878870894914371762771431462665764010129192912372490340449901", "5239065275003145843160321807696531775964858360555566589197008236687533209496"],
        ["7985258549919592662769781896447490440621354347569971700598437766156081995625", "9376351072866485300578251734844671764089160611668390200194570180225759013543"],
        ["9570976950823929161626934660575939683401710897903342799921775980893943353035", "6407880900662180043240104510114613236916437723065414158006054747177494383655"],
        ["17962366505931708682321542383646032762931774796150042922562707170594807376009", "6245130621382842925623937534683990375669631277871468906941032622563934866013"],
    ];
    
    let round_keys_y_strings = [
        ["8755297148735710088898562298102910035419345760166413737479281674630323398284", "16133435893292874812888083849160666046321318009323051176910097996974633748758"],
        ["5240474505904316858775051800099222288270827863409873986701694203345984265770", "16516377322346822856154252461095180562000423191949949242508439100972699801595"],
        ["9012679925958717565787111885188464538194947839997341443807348023221726055342", "3513323292129390671339287145562649862242777741759770715956300048086055264273"],
        ["21855834035835287540286238525800162342051591799629360593177152465113152235615", "5945179541709432313351711573896685950772105367183734375093638912196647730870"],
        ["11227229470941648605622822052481187204980748641142847464327016901091886692935", "874490282529106871250179638055108647411431264552976943414386206857408624500"],
        ["8277823808153992786803029269162651355418392229624501612473854822154276610437", "14911320361190879980016686915823914584756893340104182663424627943175208757859"],
        ["20904607884889140694334069064199005451741168419308859136555043894134683701950", "15657880601171476575713502187548665287918791967520790431542060879010363657805"],
        ["1902748146936068574869616392736208205391158973416079524055965306829204527070", "14311738005510898661766244714944477794557156116636816483240167459479765463026"],
        ["14452570815461138929654743535323908350592751448372202277464697056225242868484", "18878429879072656191963192145256996413709289475622337294803628783509021017215"],
        ["10548134661912479705005015677785100436776982856523954428067830720054853946467", "21613568037783775488400147863112554980555854603176833550688470336449256480025"],
        ["17068729307795998980462158858164249718900656779672000551618940554342475266265", "2490802518193809975066473675670874471230712567215812226164489400543194289596"],
        ["16199718037005378969178070485166950928725365516399196926532630556982133691321", "21217120779706380859547833993003263088538196273665904984368420139631145468592"],
        ["19148564379197615165212957504107910110246052442686857059768087896511716255278", "19611778548789975299387421023085714500105803761017217976092023831374602045251"],
        ["5497141763311860520411283868772341077137612389285480008601414949457218086902", "19294458970356379238521378434506704614768857764591229894917601756581488831876"],
    ];
    
    // PREPROCESSED_ROUND_KEYS for qprk values (these are the actual selector values)
    let preprocessed_round_keys_x_strings = [
        ["9875235397644879082677551174832367614794066768374461301425281161472772669364", "21858645442887666000649962444987448281406846313183347319591597416372520936186"],
        ["16782726861879113354475406688883165555010923651788393550320367281500279757725", "21716573900346641246759819543810812868751270056692513002115190154768694264248"],
        ["21061966870155710578694816572762821778055453072317217584284979102445184722013", "4549699248331629386178256801656616048620437157601772274159045439897291365892"],
        ["2305171087224367456066076264392784488317998917489739225874252016996953925819", "753423837773794493583073069146671816131950370427940330023762710613130114284"],
        ["10698642504920643104042109480794649581252624576357861152208736776268689349977", "15861688217704403960557684449792840983848058583334483774016142194529689550996"],
        ["423541189543553676504695293714968243652665955997218133034093874913252335294", "2085108831563138776622770735491169518313516126676144009754728151028847503805"],
        ["9296577236668306825777791135277154873361365777577262330385210251177355111236", "13116726794489691995858453211791287240125772343837695887860897125838790404152"],
        ["7770554041004090339324059898697839176238158454150868392168895497050756225989", "15470974097886414257515842937693682793336718238774806265675536486339112492265"],
        ["2027607608657541590621559284719638781643880713344167215688101555484801034021", "17300926706027606418491168431313029799747253325902468045764528997217592945985"],
        ["15221128912303148791281337033495862121987426916739733934930599252800117682974", "13613239587442289301259781578364994509804749661959905633737037792945804710990"],
        ["12005763401209949733243841400657315836955319000163295072425874164656732641981", "22705376494938444016386495814792101413510899935471717349511104697970912694"],
        ["14955552946683837036542615852864240457364059063741425926890312579929291591324", "288970556961158641815624492891121300982160056362935242882312884201067196942"],
        ["1581177506494232880889471258607800990368585397168021361552725068313737403708", "17981970841153864433894117072118866447373490474241757839845618545859583545511"],
        ["13826749781571877127981823420956728751726595679537198790271191951241691704832", "21456423298419106344829058807152140321464760141466308742744827391066120856237"],
    ];
    
    let preprocessed_round_keys_y_strings = [
        ["13004335645468876947782817511996516830557692388848756239167689579223703209154", "11864075285365324632501660503932294097119662259150439783414276164786389548361"],
        ["7862495485034485030006053329979953690634378679977822019470434513025641948468", "21778305909519758427732388472275509688429815430677336887808798633722753860845"],
        ["12931102024200069317238425826876622077088120606615813415940805446744126635881", "7837661096836606004314569173269958689435480650877562806125674115879275837181"],
        ["14988274660376568290931678743130590897577302840578069596030418254228064426148", "14818345905108638164688641616372585080538194792946547345972306256783653004291"],
        ["11966397199239721279456793945370572038235535122896503364930899557716957223959", "2853352834541474475776137785488700155363788984994460875907827022572233875584"],
        ["6473747423912923573021858532418794714202395821695920085715793777854113577052", "14603107593725024233314048684876188310197903996220843563409821502003190608529"],
        ["10018141451544555380964804958767235988622088919781088363105734833991047400353", "83445762062875740982996603123888928543771735703494814377210678846969285492"],
        ["4853894954678028326595990416033041454601372518726024075995342652050367914374", "13529950793291157200862531998635554831775405064348392294420225414209107516565"],
        ["2807960038839395770936423062783538297061734914581689261511199436908401205594", "2959287061458222329954767340179788517820610776269329086252152135975612854535"],
        ["1011199386146110957860470152252409466117369100436101125582703221610204956433", "11916226082408980147601015423483352131731691070197152337036758512414772646884"],
        ["6143620485513326860817743193060169274247928932037486340946124795304534640217", "9249411266687228683218385140647077689008430999582930158150164394401064685612"],
        ["3865024776110368315374386772707941373393630458661571912715432285796031519218", "1124707246745155402135444593036779382685949620543378002907953793036433309720"],
        ["3747281796037953947554825537973345299481414684769676458997083724683939123632", "516368516371014501734378213574699667472834788768584825363152893957289265859"],
        ["8415215912404504262033404854405294287543393294861880020399730040978826989992", "10041866203038674311709434184968252713427481340634280330144689403763670911641"],
    ];
    
    // Convert string constants to F
    let round_keys_x: Vec<[F; 2]> = round_keys_x_strings
        .iter()
        .map(|round_key| [
            F::from_str_vartime(round_key[0]).unwrap_or(F::ZERO),
            F::from_str_vartime(round_key[1]).unwrap_or(F::ZERO),
        ])
        .collect();
    
    let round_keys_y: Vec<[F; 2]> = round_keys_y_strings
        .iter()
        .map(|round_key| [
            F::from_str_vartime(round_key[0]).unwrap_or(F::ZERO),
            F::from_str_vartime(round_key[1]).unwrap_or(F::ZERO),
        ])
        .collect();
    
    // Convert preprocessed round keys for qprk values
    let preprocessed_round_keys_x: Vec<[F; 2]> = preprocessed_round_keys_x_strings
        .iter()
        .map(|round_key| [
            F::from_str_vartime(round_key[0]).unwrap_or(F::ZERO),
            F::from_str_vartime(round_key[1]).unwrap_or(F::ZERO),
        ])
        .collect();
    
    let preprocessed_round_keys_y: Vec<[F; 2]> = preprocessed_round_keys_y_strings
        .iter()
        .map(|round_key| [
            F::from_str_vartime(round_key[0]).unwrap_or(F::ZERO),
            F::from_str_vartime(round_key[1]).unwrap_or(F::ZERO),
        ])
        .collect();

    // println!("preprocessed_round_keys_x: {:?}", preprocessed_round_keys_x);
    // println!("preprocessed_round_keys_y: {:?}", preprocessed_round_keys_y);
    
    // Alpha inverse values for S-box
    let alpha_inv = [14981214993055009997u64, 6006880321387387405u64, 10624953561019755799u64, 2789598613442376532u64];
    
    // Create initial state (2 field elements each for x and y)
    // let mut current_x = [F::random(&mut witness_rng), F::random(&mut witness_rng)];
    // let mut current_y = [F::random(&mut witness_rng), F::ZERO]; // Salt is 0
    let mut current_x = [F::from(1u64), F::from(2u64)];
    let mut current_y = [F::from(3u64), F::ZERO];
    
    let actual_gates = GATES_NEEDED.min(usable_indices.len());
    
    // Input gate
    if !usable_indices.is_empty() {
        let idx = usable_indices[0];
        w1_values[idx] = current_x[0];
        w2_values[idx] = current_x[1];
        w3_values[idx] = current_y[0];
        w4_values[idx] = current_y[1];
        wo_values[idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1]; // Jive output
        
        // Simple constraint for input gate
        q1[idx] = F::ONE;
        q2[idx] = F::ONE;
        q3[idx] = F::ONE;
        q4[idx] = F::ONE;
        qo[idx] = F::ONE;
    }
    
    // Anemoi rounds - simulate the actual Anemoi permutation following mod.rs
    for round in 0..NUM_ROUNDS {
        if round + 1 >= actual_gates {
            break;
        }
        
        let idx = usable_indices[round + 1];
        
        // Store pre-round values
        w1_values[idx] = current_x[0];
        w2_values[idx] = current_x[1];
        w3_values[idx] = current_y[0];
        w4_values[idx] = current_y[1];
        
        // Apply one round of Anemoi permutation following mod.rs implementation
        // Step 1: Add round constants
        current_x[0] += round_keys_x[round][0];
        current_x[1] += round_keys_x[round][1];
        current_y[0] += round_keys_y[round][0];
        current_y[1] += round_keys_y[round][1];
        
        // Step 2: MDS matrix application
        let temp_x0 = current_x[0] + generator * current_x[1];
        let temp_x1 = generator * current_x[0] + generator_square_plus_one * current_x[1];
        current_x[0] = temp_x0;
        current_x[1] = temp_x1;
        
        // MDS to y with word permutation (y[1], y[0] -> y[0], y[1])
        let temp_y0 = current_y[1] + generator * current_y[0];
        let temp_y1 = generator * current_y[1] + generator_square_plus_one * current_y[0];
        current_y[0] = temp_y0;
        current_y[1] = temp_y1;
        
        // Step 3: PHT transformation (if enabled)
        for i in 0..2 {
            current_y[i] += current_x[i];
            current_x[i] += current_y[i];
        }
        
        // Step 4: S-box application
        for i in 0..2 {
            current_x[i] -= generator * current_y[i].square();
            current_y[i] = current_y[i].pow(alpha_inv); // Use alpha_inv for power
            current_x[i] += generator * current_y[i].square() + generator_inv;
        }
        
        // Set preprocessed round key selectors (use preprocessed values for qprk)
        qprk1[idx] = preprocessed_round_keys_x[round][0];
        qprk2[idx] = preprocessed_round_keys_x[round][1];
        qprk3[idx] = preprocessed_round_keys_y[round][0];
        qprk4[idx] = preprocessed_round_keys_y[round][1];
        
        // Output for this round
        wo_values[idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1];
        
        // Only add copy constraint for wo -> w4 connection as you mentioned
        if round > 0 {
            let prev_idx = usable_indices[round];
            permutation.copy((4, prev_idx), (3, idx)); // wo -> w4 only
        }
    }
    
    // Final transformations (post-rounds) following mod.rs
    // Final MDS application
    let temp_x0 = current_x[0] + generator * current_x[1];
    let temp_x1 = generator * current_x[0] + generator_square_plus_one * current_x[1];
    current_x[0] = temp_x0;
    current_x[1] = temp_x1;
    
    let temp_y0 = current_y[1] + generator * current_y[0];
    let temp_y1 = generator * current_y[1] + generator_square_plus_one * current_y[0];
    current_y[0] = temp_y0;
    current_y[1] = temp_y1;
    
    // Final PHT transformation
    for i in 0..2 {
        current_y[i] += current_x[i];
        current_x[i] += current_y[i];
    }
    
    // Final output gate
    if NUM_ROUNDS + 1 < actual_gates {
        let idx = usable_indices[NUM_ROUNDS + 1];
        w1_values[idx] = current_x[0];
        w2_values[idx] = current_x[1];
        w3_values[idx] = current_y[0];
        w4_values[idx] = current_y[1];
        wo_values[idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1]; // Final Jive output
        
        // Simple constraint for final output
        q1[idx] = F::ONE;
        q2[idx] = F::ONE;
        q3[idx] = F::ONE;
        q4[idx] = F::ONE;
        qo[idx] = F::ONE;
    }
    
    println!("w1_values: {:?}", w1_values);
    println!("w2_values: {:?}", w2_values);
    println!("w3_values: {:?}", w3_values);
    println!("w4_values: {:?}", w4_values);
    println!("wo_values: {:?}", wo_values);


    // Create circuit info
    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![], // No instances
        preprocess_polys: vec![
            q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb,
            qprk1, qprk2, qprk3, qprk4
        ],
        num_witness_polys: vec![5],
        num_challenges: vec![0],
        constraints: anemoi_hash_circuit_info::<F>(num_vars, NUM_ROUNDS).constraints,
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(7),
    };
    
    let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];
    
    (
        circuit_info,
        MockCircuit::new(vec![], witness),
    )
}





// pub fn rand_anemoi_hash_circuit<F: PrimeField, R: Rotatable + From<usize>>(
//     mut _preprocess_rng: impl RngCore,
//     mut witness_rng: impl RngCore,
// ) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
//     const NUM_ROUNDS: usize = 14; // From AnemoiJive256 specification
//     const GATES_NEEDED: usize = NUM_ROUNDS + 3; // Input + rounds + output + padding
    
//     // Calculate minimum num_vars needed
//     let num_vars = (GATES_NEEDED as f64).log2().ceil() as usize;
//     let size = 1 << num_vars;
//     let usable_indices = R::from(num_vars).usable_indices();
    
//     // Initialize polynomials
//     let mut w1_values = vec![F::ZERO; size];
//     let mut w2_values = vec![F::ZERO; size];
//     let mut w3_values = vec![F::ZERO; size];
//     let mut w4_values = vec![F::ZERO; size];
//     let mut wo_values = vec![F::ZERO; size];
    
//     // Selector polynomials (14 total)
//     let mut q1 = vec![F::ZERO; size];
//     let mut q2 = vec![F::ZERO; size];
//     let mut q3 = vec![F::ZERO; size];
//     let mut q4 = vec![F::ZERO; size];
//     let mut qo = vec![F::ZERO; size];
//     let mut qm1 = vec![F::ZERO; size];
//     let mut qm2 = vec![F::ZERO; size];
//     let mut qc = vec![F::ZERO; size];
//     let mut qecc = vec![F::ZERO; size];
//     let mut qb = vec![F::ZERO; size];
//     let mut qprk1 = vec![F::ZERO; size];
//     let mut qprk2 = vec![F::ZERO; size];
//     let mut qprk3 = vec![F::ZERO; size];
//     let mut qprk4 = vec![F::ZERO; size];
    
//     let mut permutation = Permutation::default();
    
//     // Use real Anemoi constants
//     let generator = F::from(5u64);
//     let generator_inv = F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap_or(generator.invert().unwrap());
//     let generator_square_plus_one = F::from(26u64);
    
//     // Real ROUND_KEYS_X from AnemoiJive254
//     let round_keys_x_strings = [
//         ["37", "3751828524803055471428227881618625174556947755988347881191159153764975591158"],
//         ["13352247125433170118601974521234241686699252132838635793584252509352796067497", "21001839722121566863419881512791069124083822968210421491151340238400176843969"],
//         ["8959866518978803666083663798535154543742217570455117599799616562379347639707", "21722442537234642741320951134727484119993387379465291657407115605240150584902"],
//         ["3222831896788299315979047232033900743869692917288857580060845801753443388885", "5574110054747610058729632355948568604793546392090976147435879266833412620404"],
//         ["11437915391085696126542499325791687418764799800375359697173212755436799377493", "19347108854758320361854968987183753113398822331033233961719129079198795045322"],
//         ["14725846076402186085242174266911981167870784841637418717042290211288365715997", "17733032409684964025894538244134113560864261458948810209753406163729963104066"],
//         ["3625896738440557179745980526949999799504652863693655156640745358188128872126", "16641102106808059030810525726117803887885616319153331237086309361060282564245"],
//         ["463291105983501380924034618222275689104775247665779333141206049632645736639", "9245970744804222215259369270991414441925747897718226734085751033703871913242"],
//         ["17443852951621246980363565040958781632244400021738903729528591709655537559937", "18243401795478654990110719981452738859015913555820749188627866268359980949315"],
//         ["10761214205488034344706216213805155745482379858424137060372633423069634639664", "18200337361605220875540054729693479452916227111908726624753615870884702413869"],
//         ["1555059412520168878870894914371762771431462665764010129192912372490340449901", "5239065275003145843160321807696531775964858360555566589197008236687533209496"],
//         ["7985258549919592662769781896447490440621354347569971700598437766156081995625", "9376351072866485300578251734844671764089160611668390200194570180225759013543"],
//         ["9570976950823929161626934660575939683401710897903342799921775980893943353035", "6407880900662180043240104510114613236916437723065414158006054747177494383655"],
//         ["17962366505931708682321542383646032762931774796150042922562707170594807376009", "6245130621382842925623937534683990375669631277871468906941032622563934866013"],
//     ];
    
//     // Real ROUND_KEYS_Y from AnemoiJive254
//     let round_keys_y_strings = [
//         ["8755297148735710088898562298102910035419345760166413737479281674630323398284", "16133435893292874812888083849160666046321318009323051176910097996974633748758"],
//         ["5240474505904316858775051800099222288270827863409873986701694203345984265770", "16516377322346822856154252461095180562000423191949949242508439100972699801595"],
//         ["9012679925958717565787111885188464538194947839997341443807348023221726055342", "3513323292129390671339287145562649862242777741759770715956300048086055264273"],
//         ["21855834035835287540286238525800162342051591799629360593177152465113152235615", "5945179541709432313351711573896685950772105367183734375093638912196647730870"],
//         ["11227229470941648605622822052481187204980748641142847464327016901091886692935", "874490282529106871250179638055108647411431264552976943414386206857408624500"],
//         ["8277823808153992786803029269162651355418392229624501612473854822154276610437", "14911320361190879980016686915823914584756893340104182663424627943175208757859"],
//         ["20904607884889140694334069064199005451741168419308859136555043894134683701950", "15657880601171476575713502187548665287918791967520790431542060879010363657805"],
//         ["1902748146936068574869616392736208205391158973416079524055965306829204527070", "14311738005510898661766244714944477794557156116636816483240167459479765463026"],
//         ["14452570815461138929654743535323908350592751448372202277464697056225242868484", "18878429879072656191963192145256996413709289475622337294803628783509021017215"],
//         ["10548134661912479705005015677785100436776982856523954428067830720054853946467", "21613568037783775488400147863112554980555854603176833550688470336449256480025"],
//         ["17068729307795998980462158858164249718900656779672000551618940554342475266265", "2490802518193809975066473675670874471230712567215812226164489400543194289596"],
//         ["16199718037005378969178070485166950928725365516399196926532630556982133691321", "21217120779706380859547833993003263088538196273665904984368420139631145468592"],
//         ["19148564379197615165212957504107910110246052442686857059768087896511716255278", "19611778548789975299387421023085714500105803761017217976092023831374602045251"],
//         ["5497141763311860520411283868772341077137612389285480008601414949457218086902", "19294458970356379238521378434506704614768857764591229894917601756581488831876"],
//     ];
    
//     // Convert string constants to F
//     let round_keys_x: Vec<[F; 2]> = round_keys_x_strings
//         .iter()
//         .map(|round_key| [
//             F::from_str_vartime(round_key[0]).unwrap_or(F::ZERO),
//             F::from_str_vartime(round_key[1]).unwrap_or(F::ZERO),
//         ])
//         .collect();
    
//     let round_keys_y: Vec<[F; 2]> = round_keys_y_strings
//         .iter()
//         .map(|round_key| [
//             F::from_str_vartime(round_key[0]).unwrap_or(F::ZERO),
//             F::from_str_vartime(round_key[1]).unwrap_or(F::ZERO),
//         ])
//         .collect();
    
//     // Create initial state (2 field elements each for x and y)
//     let mut current_x = [F::random(&mut witness_rng), F::random(&mut witness_rng)];
//     let mut current_y = [F::random(&mut witness_rng), F::ZERO]; // Salt is 0

//     println!("current_x: {:?}", current_x);
//     println!("current_y: {:?}", current_y);
    
//     let actual_gates = GATES_NEEDED.min(usable_indices.len());
    
//     // Input gate
//     if !usable_indices.is_empty() {
//         let idx = usable_indices[0];
//         w1_values[idx] = current_x[0];
//         w2_values[idx] = current_x[1];
//         w3_values[idx] = current_y[0];
//         w4_values[idx] = current_y[1];
//         wo_values[idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1]; // Jive output
        
//         // Simple constraint for input gate
//         q1[idx] = F::ONE;
//         q2[idx] = F::ONE;
//         q3[idx] = F::ONE;
//         q4[idx] = F::ONE;
//         qo[idx] = F::ONE;
//     }
    
//     // Anemoi rounds - simulate the actual Anemoi permutation
//     for round in 0..NUM_ROUNDS {
//         if round + 1 >= actual_gates {
//             break;
//         }
        
//         let idx = usable_indices[round + 1];
        
//         // Get round constants (already converted to F)
//         let round_x = if round < round_keys_x.len() {
//             round_keys_x[round]
//         } else {
//             [F::ZERO, F::ZERO]
//         };
//         let round_y = if round < round_keys_y.len() {
//             round_keys_y[round]
//         } else {
//             [F::ZERO, F::ZERO]
//         };
        
//         // Store pre-round values
//         w1_values[idx] = current_x[0];
//         w2_values[idx] = current_x[1];
//         w3_values[idx] = current_y[0];
//         w4_values[idx] = current_y[1];
        
//         // Apply one round of Anemoi permutation
//         // Add round constants
//         current_x[0] += round_x[0];
//         current_x[1] += round_x[1];
//         current_y[0] += round_y[0];
//         current_y[1] += round_y[1];
        
//         // Apply MDS matrix using converted constants
//         let temp_x0 = current_x[0] + generator * current_x[1];
//         let temp_x1 = generator * current_x[0] + generator_square_plus_one * current_x[1];
//         current_x[0] = temp_x0;
//         current_x[1] = temp_x1;
        
//         // Apply MDS to y with word permutation
//         let temp_y0 = current_y[1] + generator * current_y[0];
//         let temp_y1 = generator * current_y[1] + generator_square_plus_one * current_y[0];
//         current_y[0] = temp_y0;
//         current_y[1] = temp_y1;
        
//         // PHT transformation (if used)
//         if AnemoiJive256::USE_PHT {
//             for i in 0..2 {
//                 current_y[i] += current_x[i];
//                 current_x[i] += current_y[i];
//             }
//         }
        
//         // Apply S-box
//         for i in 0..2 {
//             current_x[i] -= generator * current_y[i].square();
//             current_y[i] = current_y[i].pow([AnemoiJive256::ALPHA as u64]);
//             current_x[i] += generator * current_y[i].square() + generator_inv;
//         }
        
//         // Set preprocessed round key selectors (use converted values)
//         qprk1[idx] = round_x[0];
//         qprk2[idx] = round_x[1];
//         qprk3[idx] = round_y[0];
//         qprk4[idx] = round_y[1];
        
//         // Output for this round
//         wo_values[idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1];
        
//         // Only add copy constraint for wo -> w4 connection as you mentioned
//         if round > 0 {
//             let prev_idx = usable_indices[round];
//             permutation.copy((4, prev_idx), (3, idx)); // wo -> w4 only
//         }
//     }
    
//     // Final MDS application (post-rounds)
//     let temp_x0 = current_x[0] + generator * current_x[1];
//     let temp_x1 = generator * current_x[0] + generator_square_plus_one * current_x[1];
//     current_x[0] = temp_x0;
//     current_x[1] = temp_x1;
    
//     let temp_y0 = current_y[1] + generator * current_y[0];
//     let temp_y1 = generator * current_y[1] + generator_square_plus_one * current_y[0];
//     current_y[0] = temp_y0;
//     current_y[1] = temp_y1;
    
//     if AnemoiJive256::USE_PHT {
//         for i in 0..2 {
//             current_y[i] += current_x[i];
//             current_x[i] += current_y[i];
//         }
//     }
    
//     // Final output gate
//     if NUM_ROUNDS + 1 < actual_gates {
//         let idx = usable_indices[NUM_ROUNDS + 1];
//         w1_values[idx] = current_x[0];
//         w2_values[idx] = current_x[1];
//         w3_values[idx] = current_y[0];
//         w4_values[idx] = current_y[1];
//         wo_values[idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1]; // Final Jive output
        
//         // Simple constraint for final output
//         q1[idx] = F::ONE;
//         q2[idx] = F::ONE;
//         q3[idx] = F::ONE;
//         q4[idx] = F::ONE;
//         qo[idx] = F::ONE;
//     }
    
//     println!("wo_values: {:?}", w1_values);

    
//     // Create circuit info
//     let circuit_info = PlonkishCircuitInfo {
//         k: num_vars,
//         num_instances: vec![],
//         preprocess_polys: vec![
//             q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb,
//             qprk1, qprk2, qprk3, qprk4
//         ],
//         num_witness_polys: vec![5],
//         num_challenges: vec![0],
//         constraints: anemoi_hash_circuit_info::<F>(num_vars, NUM_ROUNDS).constraints,
//         lookups: Vec::new(),
//         permutations: permutation.into_cycles(),
//         max_degree: Some(7),
//     };
    
//     let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];
    
//     (
//         circuit_info,
//         MockCircuit::new(vec![], witness),
//     )
// }

// /// Create a circuit that computes Anemoi Jive compression (3-to-1)
// pub fn anemoi_jive_compression_circuit<F: PrimeField>(
//     num_vars: usize,
// ) -> PlonkishCircuitInfo<F> {
//     const NUM_ROUNDS: usize = 14; // From the paper
    
//     // Create base Anemoi circuit
//     let mut circuit_info = anemoi_hash_circuit_info::<F>(num_vars, NUM_ROUNDS);
    
//     // Add constraint for final sum computation (Jive mode)
//     // The output is sum of input and output of permutation
//     let output_sum_gate = NUM_ROUNDS * 16; // After all round gates
    
//     // Additional constraints for Jive compression can be added here
    
//     circuit_info
// }


pub fn vanilla_plonk_circuit_info<F: PrimeField>(
    num_vars: usize,
    num_instances: usize,
    preprocess_polys: [Vec<F>; 5],
    permutations: Vec<Vec<(usize, usize)>>,
) -> PlonkishCircuitInfo<F> {
    let [pi, q_l, q_r, q_m, q_o, q_c, w_l, w_r, w_o]: &[Expression<F>; 9] =
        &array::from_fn(|poly| Query::new(poly, Rotation::cur())).map(Expression::Polynomial);
    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![num_instances],
        preprocess_polys: preprocess_polys.to_vec(),
        num_witness_polys: vec![3],
        num_challenges: vec![0],
        constraints: vec![q_l * w_l + q_r * w_r + q_m * w_l * w_r + q_o * w_o + q_c + pi],
        lookups: Vec::new(),
        permutations,
        max_degree: Some(4),
    }
}

pub fn vanilla_plonk_expression<F: PrimeField>(num_vars: usize) -> Expression<F> {
    let circuit_info = vanilla_plonk_circuit_info(
        num_vars,
        0,
        Default::default(),
        vec![vec![(6, 1)], vec![(7, 1)], vec![(8, 1)]],
    );
    let (num_permutation_z_polys, expression) = compose(&circuit_info);
    assert_eq!(num_permutation_z_polys, 1);
    expression
}

pub fn vanilla_plonk_w_lookup_circuit_info<F: PrimeField>(
    num_vars: usize,
    num_instances: usize,
    preprocess_polys: [Vec<F>; 9],
    permutations: Vec<Vec<(usize, usize)>>,
) -> PlonkishCircuitInfo<F> {
    let [pi, q_l, q_r, q_m, q_o, q_c, q_lookup, t_l, t_r, t_o, w_l, w_r, w_o] =
        &array::from_fn(|poly| Query::new(poly, Rotation::cur())).map(Expression::Polynomial);
    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![num_instances],
        preprocess_polys: preprocess_polys.to_vec(),
        num_witness_polys: vec![3],
        num_challenges: vec![0],
        constraints: vec![q_l * w_l + q_r * w_r + q_m * w_l * w_r + q_o * w_o + q_c + pi],
        lookups: vec![vec![
            (q_lookup * w_l, t_l.clone()),
            (q_lookup * w_r, t_r.clone()),
            (q_lookup * w_o, t_o.clone()),
        ]],
        permutations,
        max_degree: Some(4),
    }
}

pub fn vanilla_plonk_w_lookup_expression<F: PrimeField>(num_vars: usize) -> Expression<F> {
    let circuit_info = vanilla_plonk_w_lookup_circuit_info(
        num_vars,
        0,
        Default::default(),
        vec![vec![(10, 1)], vec![(11, 1)], vec![(12, 1)]],
    );
    let (num_permutation_z_polys, expression) = compose(&circuit_info);
    assert_eq!(num_permutation_z_polys, 1);
    expression
}

pub fn rand_vanilla_plonk_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    num_vars: usize,
    mut preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    let size = 1 << num_vars;
    let mut polys = [(); 9].map(|_| vec![F::ZERO; size]);

    let instances = rand_vec(num_vars, &mut witness_rng);
    polys[0] = mem::take(&mut instance_polys::<_, R>(num_vars, [&instances])[0]).into_evals();

    let mut permutation = Permutation::default();
    for poly in [6, 7, 8] {
        permutation.copy((poly, 1), (poly, 1));
    }
    for idx in 0..size - 1 {
        let [w_l, w_r] = if preprocess_rng.next_u32().is_even() && idx > 1 {
            let [l_copy_idx, r_copy_idx] = [(); 2].map(|_| {
                (
                    rand_idx(6..9, &mut preprocess_rng),
                    rand_idx(1..idx, &mut preprocess_rng),
                )
            });
            permutation.copy(l_copy_idx, (6, idx));
            permutation.copy(r_copy_idx, (7, idx));
            [
                polys[l_copy_idx.0][l_copy_idx.1],
                polys[r_copy_idx.0][r_copy_idx.1],
            ]
        } else {
            rand_array(&mut witness_rng)
        };
        let q_c = F::random(&mut preprocess_rng);
        let values = if preprocess_rng.next_u32().is_even() {
            vec![
                (1, F::ONE),
                (2, F::ONE),
                (4, -F::ONE),
                (5, q_c),
                (6, w_l),
                (7, w_r),
                (8, w_l + w_r + q_c + polys[0][idx]),
            ]
        } else {
            vec![
                (3, F::ONE),
                (4, -F::ONE),
                (5, q_c),
                (6, w_l),
                (7, w_r),
                (8, w_l * w_r + q_c + polys[0][idx]),
            ]
        };
        for (poly, value) in values {
            polys[poly][idx] = value;
        }
    }

    let [_, q_l, q_r, q_m, q_o, q_c, w_l, w_r, w_o] = polys;
    let circuit_info = vanilla_plonk_circuit_info(
        num_vars,
        instances.len(),
        [q_l, q_r, q_m, q_o, q_c],
        permutation.into_cycles(),
    );
    (
        circuit_info,
        MockCircuit::new(vec![instances], vec![w_l, w_r, w_o]),
    )
}

pub fn rand_vanilla_plonk_assignment<F: PrimeField, R: Rotatable + From<usize>>(
    num_vars: usize,
    mut preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (Vec<MultilinearPolynomial<F>>, Vec<F>) {
    let (polys, permutations) = {
        let (circuit_info, circuit) =
            rand_vanilla_plonk_circuit::<_, R>(num_vars, &mut preprocess_rng, &mut witness_rng);
        let witness = circuit.synthesize(0, &[]).unwrap();
        let polys = chain![
            instance_polys::<_, R>(num_vars, circuit.instances()),
            chain![circuit_info.preprocess_polys, witness].map(MultilinearPolynomial::new),
        ]
        .collect_vec();
        (polys, circuit_info.permutations)
    };
    let challenges: [_; 3] = rand_array(&mut witness_rng);
    let [beta, gamma, _] = challenges;

    let permutation_polys = permutation_polys(num_vars, &[6, 7, 8], &permutations);
    let permutation_z_polys = permutation_z_polys::<_, R>(
        1,
        &[6, 7, 8]
            .into_iter()
            .zip(permutation_polys.iter().cloned())
            .collect_vec(),
        &polys.iter().collect_vec(),
        &beta,
        &gamma,
    );

    (
        chain![polys, permutation_polys, permutation_z_polys].collect_vec(),
        challenges.to_vec(),
    )
}

pub fn rand_vanilla_plonk_w_lookup_circuit<F: PrimeField + Hash, R: Rotatable + From<usize>>(
    num_vars: usize,
    mut preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    let size = 1 << num_vars;
    let mut polys = [(); 13].map(|_| vec![F::ZERO; size]);

    let [t_l, t_r, t_o] = [(); 3].map(|_| {
        chain![
            [F::ZERO; 2],
            iter::repeat_with(|| F::random(&mut preprocess_rng)),
        ]
        .take(size)
        .collect_vec()
    });
    polys[7] = t_l;
    polys[8] = t_r;
    polys[9] = t_o;

    let instances = rand_vec(num_vars, &mut witness_rng);
    polys[0] = instance_polys::<_, R>(num_vars, [&instances])[0]
        .evals()
        .to_vec();
    let instance_rows = R::from(num_vars)
        .usable_indices()
        .into_iter()
        .take(num_vars + 1)
        .collect::<HashSet<_>>();

    let mut permutation = Permutation::default();
    for poly in [10, 11, 12] {
        permutation.copy((poly, 1), (poly, 1));
    }
    for idx in 0..size - 1 {
        let use_copy = preprocess_rng.next_u32().is_even() && idx > 1;
        let [w_l, w_r] = if use_copy {
            let [l_copy_idx, r_copy_idx] = [(); 2].map(|_| {
                (
                    rand_idx(10..13, &mut preprocess_rng),
                    rand_idx(1..idx, &mut preprocess_rng),
                )
            });
            permutation.copy(l_copy_idx, (10, idx));
            permutation.copy(r_copy_idx, (11, idx));
            [
                polys[l_copy_idx.0][l_copy_idx.1],
                polys[r_copy_idx.0][r_copy_idx.1],
            ]
        } else {
            rand_array(&mut witness_rng)
        };
        let q_c = F::random(&mut preprocess_rng);
        let values = match (
            use_copy || instance_rows.contains(&idx),
            preprocess_rng.next_u32().is_even(),
        ) {
            (true, true) => {
                vec![
                    (1, F::ONE),
                    (2, F::ONE),
                    (4, -F::ONE),
                    (5, q_c),
                    (10, w_l),
                    (11, w_r),
                    (12, w_l + w_r + q_c + polys[0][idx]),
                ]
            }
            (true, false) => {
                vec![
                    (3, F::ONE),
                    (4, -F::ONE),
                    (5, q_c),
                    (10, w_l),
                    (11, w_r),
                    (12, w_l * w_r + q_c + polys[0][idx]),
                ]
            }
            (false, _) => {
                let idx = rand_idx(1..size, &mut witness_rng);
                vec![
                    (6, F::ONE),
                    (10, polys[7][idx]),
                    (11, polys[8][idx]),
                    (12, polys[9][idx]),
                ]
            }
        };
        for (poly, value) in values {
            polys[poly][idx] = value;
        }
    }

    let [_, q_l, q_r, q_m, q_o, q_c, q_lookup, t_l, t_r, t_o, w_l, w_r, w_o] = polys;
    let circuit_info = vanilla_plonk_w_lookup_circuit_info(
        num_vars,
        instances.len(),
        [q_l, q_r, q_m, q_o, q_c, q_lookup, t_l, t_r, t_o],
        permutation.into_cycles(),
    );
    (
        circuit_info,
        MockCircuit::new(vec![instances], vec![w_l, w_r, w_o]),
    )
}

pub fn rand_vanilla_plonk_w_lookup_assignment<F: PrimeField + Hash, R: Rotatable + From<usize>>(
    num_vars: usize,
    mut preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (Vec<MultilinearPolynomial<F>>, Vec<F>) {
    let (polys, permutations) = {
        let (circuit_info, circuit) = rand_vanilla_plonk_w_lookup_circuit::<_, R>(
            num_vars,
            &mut preprocess_rng,
            &mut witness_rng,
        );
        let witness = circuit.synthesize(0, &[]).unwrap();
        let polys = chain![
            instance_polys::<_, R>(num_vars, circuit.instances()),
            chain![circuit_info.preprocess_polys, witness].map(MultilinearPolynomial::new),
        ]
        .collect_vec();
        (polys, circuit_info.permutations)
    };
    let challenges: [_; 3] = rand_array(&mut witness_rng);
    let [beta, gamma, _] = challenges;

    let (lookup_compressed_polys, lookup_m_polys) = {
        let PlonkishCircuitInfo { lookups, .. } =
            vanilla_plonk_w_lookup_circuit_info(0, 0, Default::default(), Vec::new());
        let polys = polys.iter().collect_vec();
        let betas = powers(beta).take(3).collect_vec();
        let lookup_compressed_polys =
            lookup_compressed_polys::<_, R>(&lookups, &polys, &[], &betas);
        let lookup_m_polys = lookup_m_polys(&lookup_compressed_polys).unwrap();
        (lookup_compressed_polys, lookup_m_polys)
    };
    let lookup_h_polys = lookup_h_polys(&lookup_compressed_polys, &lookup_m_polys, &gamma);

    let permutation_polys = permutation_polys(num_vars, &[10, 11, 12], &permutations);
    let permutation_z_polys = permutation_z_polys::<_, R>(
        1,
        &[10, 11, 12]
            .into_iter()
            .zip(permutation_polys.iter().cloned())
            .collect_vec(),
        &polys.iter().collect_vec(),
        &beta,
        &gamma,
    );

    (
        chain![
            polys,
            permutation_polys,
            lookup_m_polys,
            lookup_h_polys,
            permutation_z_polys,
        ]
        .collect_vec(),
        challenges.to_vec(),
    )
}

#[derive(Default)]
pub struct Permutation {
    cycles: Vec<HashSet<(usize, usize)>>,
    cycle_idx: HashMap<(usize, usize), usize>,
}

impl Permutation {
    pub fn copy(&mut self, lhs: (usize, usize), rhs: (usize, usize)) {
        match self.cycle_idx.get(&lhs).copied() {
            Some(cycle_idx) => {
                self.cycles[cycle_idx].insert(rhs);
                self.cycle_idx.insert(rhs, cycle_idx);
            }
            None => {
                let cycle_idx = self.cycles.len();
                self.cycles.push(HashSet::from_iter([lhs, rhs]));
                for cell in [lhs, rhs] {
                    self.cycle_idx.insert(cell, cycle_idx);
                }
            }
        };
    }

    pub fn into_cycles(self) -> Vec<Vec<(usize, usize)>> {
        self.cycles
            .into_iter()
            .map(|cycle| cycle.into_iter().sorted().collect_vec())
            .collect()
    }
}


/// Helper trait extension for Expression to add pow method
trait ExpressionExt<F: Field> {
    fn pow(self, n: u32) -> Expression<F>;
}

impl<F: Field> ExpressionExt<F> for Expression<F> {
    fn pow(self, n: u32) -> Expression<F> {
        match n {
            0 => Expression::one(),
            1 => self,
            _ => {
                let half = self.clone().pow(n / 2);
                let half_squared = half.clone() * half;
                if n % 2 == 0 {
                    half_squared
                } else {
                    half_squared * self
                }
            }
        }
    }
}

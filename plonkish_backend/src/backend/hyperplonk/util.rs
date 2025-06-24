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


/// 简单的索引测试电路：验证 w_next = w_current + 1
pub fn rand_index_test_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    num_vars: usize,
    mut _preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    let size = 1 << num_vars;
    let usable_indices = R::from(num_vars).usable_indices();
    
    // 初始化多项式
    let mut q_seq = vec![F::ZERO; size];  // 选择器
    let mut w_values = vec![F::ZERO; size];  // witness
    
    let mut permutation = Permutation::default();
    permutation.copy((1, usable_indices[0]), (1, usable_indices[0])); // 自引用
    
    // 在usable_indices上设置递增序列
    for i in 0..std::cmp::min(usable_indices.len().saturating_sub(1), 8) {
        let current_idx = usable_indices[i];
        let next_idx = usable_indices[i + 1];
        
        // 设置值：0, 1, 2, 3, ...
        w_values[current_idx] = F::from(i as u64);
        w_values[next_idx] = F::from((i + 1) as u64);
        
        // 激活约束
        if i != 0 {
            q_seq[current_idx] = F::ONE;
        }
    }
    

    // 约束：q_seq * (w_next - w_current - 1) = 0
    let constraint: Expression<F> = Expression::Polynomial::<F>(Query::new(0, Rotation::cur())) *
        (Expression::Polynomial::<F>(Query::new(1, Rotation::cur())) -
         Expression::Polynomial::<F>(Query::new(1, Rotation::prev())) -
         Expression::one());
    
    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![],
        preprocess_polys: vec![q_seq],
        num_witness_polys: vec![1],
        num_challenges: vec![0],
        constraints: vec![constraint],
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(2),
    };
    
    (circuit_info, MockCircuit::new(vec![], vec![w_values]))
}




pub fn rand_jive_crh_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    input1: F,
    input2: F, 
    input3: F,
    padding_constant: F, // For domain separation
    mut _preprocess_rng: impl RngCore,
    mut _witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    const NUM_ROUNDS: usize = 14;
    const TOTAL_GATES: usize = NUM_ROUNDS + 3;
    
    let num_vars = (TOTAL_GATES as f64).log2().ceil() as usize;
    let size = 1 << num_vars;
    let usable_indices = R::from(num_vars).usable_indices();
    
    if usable_indices.len() < TOTAL_GATES {
        panic!("Not enough usable indices: {} < {}", usable_indices.len(), TOTAL_GATES);
    }
    
    // Initialize witness polynomials
    let mut w1_values = vec![F::ZERO; size];
    let mut w2_values = vec![F::ZERO; size];
    let mut w3_values = vec![F::ZERO; size];
    let mut w4_values = vec![F::ZERO; size];
    let mut wo_values = vec![F::ZERO; size];
    
    // Initialize selector polynomials
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
    let generator_inv = F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap_or(F::ZERO);
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
    
    // Alpha inverse values for S-box
    let alpha_inv = [14981214993055009997u64, 6006880321387387405u64, 10624953561019755799u64, 2789598613442376532u64];
    
    // Set initial state: x[1]=input1, x[2]=input2, y[1]=input3, y[2]=padding_constant
    let mut current_x = [input1, input2];
    let mut current_y = [input3, padding_constant];
    let sum_before_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    
    // Input gate (Gate 0) - enforce salt constraint
    let input_idx = usable_indices[0];
    w1_values[input_idx] = current_x[0];
    w2_values[input_idx] = current_x[1];
    w3_values[input_idx] = current_y[0];
    w4_values[input_idx] = current_y[1];
    wo_values[input_idx] = sum_before_perm; // Not used in input gate
    
    
    // Enforce w4 = padding_constant
    q1[input_idx] = F::ONE;
    q2[input_idx] = F::ONE;
    q3[input_idx] = F::ONE;
    q4[input_idx] = F::ONE;
    qo[input_idx] = F::ONE;
    // qc[input_idx] = -padding_constant;

    for poly_idx in 14..19 {  // w1, w2, w3, w4, wo
        permutation.copy((poly_idx, input_idx), (poly_idx, input_idx));
    }
    
    // Anemoi rounds (Gates 1-14)
    for round in 0..NUM_ROUNDS {
        let gate_idx = usable_indices[round + 1];
        
        // Set preprocessed round keys
        qprk1[gate_idx] = preprocessed_round_keys_x[round][0];
        qprk2[gate_idx] = preprocessed_round_keys_x[round][1];
        qprk3[gate_idx] = preprocessed_round_keys_y[round][0];
        qprk4[gate_idx] = preprocessed_round_keys_y[round][1];
        
        // Store pre-round state
        w1_values[gate_idx] = current_x[0];
        w2_values[gate_idx] = current_x[1];
        w3_values[gate_idx] = current_y[0];
        w4_values[gate_idx] = current_y[1];
        
        // Apply Anemoi round (reuse your existing logic)
        current_x[0] += round_keys_x[round][0];
        current_x[1] += round_keys_x[round][1];
        current_y[0] += round_keys_y[round][0];
        current_y[1] += round_keys_y[round][1];
        
        // MDS matrix
        let temp_x0 = current_x[0] + generator * current_x[1];
        let temp_x1 = generator * current_x[0] + generator_square_plus_one * current_x[1];
        current_x[0] = temp_x0;
        current_x[1] = temp_x1;
        
        let temp_y0 = current_y[1] + generator * current_y[0];
        let temp_y1 = generator * current_y[1] + generator_square_plus_one * current_y[0];
        current_y[0] = temp_y0;
        current_y[1] = temp_y1;
        
        // PHT transformation
        for i in 0..2 {
            current_y[i] += current_x[i];
            current_x[i] += current_y[i];
        }
        
        // S-box
        // let alpha_inv = [14981214993055009997u64, 6006880321387387405u64, 10624953561019755799u64, 2789598613442376532u64];
        for i in 0..2 {
            current_x[i] -= generator * current_y[i].square();
            current_y[i] -= current_x[i].pow(alpha_inv); // Use correct alpha_inv
            current_x[i] += generator * current_y[i].square() + generator_inv;
        }
    }
        // Set up next round values for the last round to enable Rotation::next() constraints
    let final_idx = usable_indices[NUM_ROUNDS + 1];
    w1_values[final_idx] = current_x[0];
    w2_values[final_idx] = current_x[1];
    w3_values[final_idx] = current_y[0];
    w4_values[final_idx] = current_y[1];

    // final state constraint
    // TODO: 这里为了安全性最好做一下拆分，目前是PoC，所以直接写在一起了
    q1[final_idx] = F::from(3u64) * generator + F::from(3u64);
    q2[final_idx] = F::from(3u64) * (generator * generator + generator + F::ONE);
    q3[final_idx] = F::from(2u64) * (generator * generator + generator + F::ONE);
    q4[final_idx] = F::from(2u64) * generator + F::from(2u64);
    qo[final_idx] = F::ONE;

    
    // Final MDS and PHT
    let temp_x0 = current_x[0] + generator * current_x[1];
    let temp_x1 = generator * current_x[0] + generator_square_plus_one * current_x[1];
    current_x[0] = temp_x0;
    current_x[1] = temp_x1;
    
    let temp_y0 = current_y[1] + generator * current_y[0];
    let temp_y1 = generator * current_y[1] + generator_square_plus_one * current_y[0];
    current_y[0] = temp_y0;
    current_y[1] = temp_y1;
    
    for i in 0..2 {
        current_y[i] += current_x[i];
        current_x[i] += current_y[i];
    }
    
    let sum_after_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    wo_values[final_idx] = sum_after_perm;
    
    // // Sum computation gate (Gate NUM_ROUNDS+1)
    // let sum_gate_idx = usable_indices[NUM_ROUNDS + 2];
    // w1_values[sum_gate_idx] = current_x[0];
    // w2_values[sum_gate_idx] = current_x[1];
    // w3_values[sum_gate_idx] = current_y[0];
    // w4_values[sum_gate_idx] = current_y[1];
    // wo_values[sum_gate_idx] = sum_after_perm;
    
    // Constraint: wo = w1 + w2 + w3 + w4 (sum of permutation output)
    // q1[sum_gate_idx] = F::ONE;
    // q2[sum_gate_idx] = F::ONE;
    // q3[sum_gate_idx] = F::ONE;
    // q4[sum_gate_idx] = F::ONE;
    // qo[sum_gate_idx] = F::ONE;
    
    // Final output gate (Gate NUM_ROUNDS+2) - Jive output
    let output_gate_idx = usable_indices[NUM_ROUNDS + 2];
    let jive_output = sum_before_perm + sum_after_perm;
    
    w1_values[output_gate_idx] = sum_before_perm;
    w2_values[output_gate_idx] = sum_after_perm;
    w3_values[output_gate_idx] = F::ZERO;
    w4_values[output_gate_idx] = F::ZERO;
    wo_values[output_gate_idx] = jive_output;
    permutation.copy((18, final_idx), (15, output_gate_idx));
    permutation.copy((14, output_gate_idx), (18, usable_indices[NUM_ROUNDS + 2 - 16]));
    
    // Constraint: wo = w1 + w2 (final Jive output)
    q1[output_gate_idx] = F::ONE;
    q2[output_gate_idx] = F::ONE;
    qo[output_gate_idx] = F::ONE;
    
    // Set up copy constraints for round transitions
    for round in 1..=NUM_ROUNDS {
        let current_idx = usable_indices[round];
        let next_idx = usable_indices[round + 1];
        
        wo_values[current_idx] = w4_values[next_idx];
        permutation.copy((18, current_idx), (17, next_idx)); // wo -> w4
    }
    
    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![],
        preprocess_polys: vec![
            q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb,
            qprk1, qprk2, qprk3, qprk4
        ],
        num_witness_polys: vec![5],
        num_challenges: vec![0],
        constraints: anemoi_hash_circuit_info::<F>(num_vars).constraints,
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(7),
    };
    
    let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];
    
    (circuit_info, MockCircuit::new(vec![], witness))
}



/// Create a PlonkishCircuitInfo for Anemoi hash with TurboPlonk constraints
/// Based on the paper "An efficient verifiable state for zk-EVM and beyond from the Anemoi hash function"
pub fn anemoi_hash_circuit_info_original<F: PrimeField>(
    num_vars: usize,
    // _num_rounds: usize,
) -> PlonkishCircuitInfo<F> {

    //14 selector polynomials for the Anemoi constraints
    // Create expressions for selector polynomials
    let [q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4] = 
        &array::from_fn(|i| Query::new(i, Rotation::cur())).map(Expression::Polynomial);

    // We need 5 wires as described in the paper: w1, w2, w3, w4, wo
    // Create expressions for wire polynomials
    let [w1, w2, w3, w4, wo] = &array::from_fn(|i| {
        Query::new(i + 14, Rotation::cur())
    }).map(Expression::Polynomial);

    
    
    // Create expressions for next rotation (for Anemoi constraints)
    let [w1_next, w2_next, w3_next, w4_next] = &[
        Query::new(14, Rotation::next()),
        Query::new(15, Rotation::next()),
        Query::new(16, Rotation::next()),
        Query::new(17, Rotation::next()),
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
    let g_inv = Expression::<F>::Constant(F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap()); // Delta
    
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
        (c_prime_2.clone() - w4_next).pow(5)
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
        (c_prime_2 - w4_next.clone()).pow(5)
        + g.clone() * w4_next.clone().pow(2)
        + g_inv.clone()
        - w2_next
    );

    let test_constraint =  g * g_inv.clone() - Expression::one() ;

    // Collect all constraints
    let mut constraints = vec![base_constraint];
    constraints.extend(bool_constraints);
    constraints.extend(vec![anemoi_1, anemoi_2, anemoi_3, anemoi_4]);
    // constraints.extend(vec![test_constraint]);
    
    
    // Create preprocessed polynomials (selectors) - 14 selectors
    let preprocess_polys = vec![vec![F::ZERO; 1 << num_vars]; 14];
    
    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![0],
        preprocess_polys,
        num_witness_polys: vec![5], // w1, w2, w3, w4, wo
        num_challenges: vec![0],
        constraints,
        // constraints: vec![test_constraint],
        lookups: Vec::new(),
        permutations: Vec::new(),
        max_degree: Some(7), // Due to pow(5) in Anemoi constraints
    }
}





/// Create a PlonkishCircuitInfo for Anemoi hash with TurboPlonk constraints
/// Based on the paper "An efficient verifiable state for zk-EVM and beyond from the Anemoi hash function"
pub fn anemoi_hash_circuit_info<F: PrimeField>(
    num_vars: usize,
    // _num_rounds: usize,
) -> PlonkishCircuitInfo<F> {

    //14 selector polynomials for the Anemoi constraints
    // Create expressions for selector polynomials
    let [q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb, qprk1, qprk2, qprk3, qprk4] = 
        &array::from_fn(|i| Query::new(i, Rotation::cur())).map(Expression::Polynomial);

    // We need 5 wires as described in the paper: w1, w2, w3, w4, wo
    // Create expressions for wire polynomials
    let [w1, w2, w3, w4, wo] = &array::from_fn(|i| {
        Query::new(i + 14, Rotation::cur())
    }).map(Expression::Polynomial);

    
    
    // Create expressions for next rotation (for Anemoi constraints)
    let [w1_next, w2_next, w3_next] = &[
        Query::new(14, Rotation::next()),
        Query::new(15, Rotation::next()),
        Query::new(16, Rotation::next()),
        // Query::new(17, Rotation::next()),
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
    let g_inv = Expression::<F>::Constant(F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap()); // Delta
    
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
        + g.clone() * wo.clone().pow(2)
        + g_inv.clone()
        - w2_next
    );

    let test_constraint =  g * g_inv.clone() - Expression::one() ;

    // Collect all constraints
    let mut constraints = vec![base_constraint];
    constraints.extend(bool_constraints);
    constraints.extend(vec![anemoi_1, anemoi_2, anemoi_3, anemoi_4]);
    // constraints.extend(vec![test_constraint]);
    
    
    // Create preprocessed polynomials (selectors) - 14 selectors
    let preprocess_polys = vec![vec![F::ZERO; 1 << num_vars]; 14];
    
    PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![0],
        preprocess_polys,
        num_witness_polys: vec![5], // w1, w2, w3, w4, wo
        num_challenges: vec![0],
        constraints,
        // constraints: vec![test_constraint],
        lookups: Vec::new(),
        permutations: Vec::new(),
        max_degree: Some(7), // Due to pow(5) in Anemoi constraints
    }
}

/// Generate a complete Anemoi hash circuit implementing the paper's specifications
/// Automatically determines the required num_vars based on gate count
/// Generate a complete Anemoi hash circuit implementing the paper's specifications
/// Automatically determines the required num_vars based on gate count
pub fn rand_anemoi_hash_circuit_with_flatten<F: PrimeField, R: Rotatable + From<usize>>(
    mut _preprocess_rng: impl RngCore,
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    const NUM_ROUNDS: usize = 14; // From AnemoiJive256 specification
    const GATES_NEEDED: usize = NUM_ROUNDS + 4; // Input + rounds + output + padding
    
    // Calculate minimum num_vars needed
    let num_vars = (GATES_NEEDED as f64).log2().ceil() as usize;
    let size = 1 << num_vars;
    
    // 获取可用的索引
    let usable_indices = R::from(num_vars).usable_indices();
    // println!("num_vars: {}, size: {}, usable_indices len: {}", num_vars, size, usable_indices.len());
    // println!("usable_indices: {:?}", &usable_indices[..std::cmp::min(usable_indices.len(), 20)]);
    
    if usable_indices.len() < GATES_NEEDED {
        panic!("Not enough usable indices: {} < {}", usable_indices.len(), GATES_NEEDED);
    }
    
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
    // 1. 初始化：为每个witness多项式创建初始cycle
    // for poly_idx in 14..=18 {  // w1, w2, w3, w4, wo
    //     permutation.copy((poly_idx, 1), (poly_idx, 1));  // 自引用
    // }

    // Use real Anemoi constants
    let generator = F::from(5u64);
    let generator_inv = F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap_or(F::ZERO);
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
    
    // Alpha inverse values for S-box
    let alpha_inv = [14981214993055009997u64, 6006880321387387405u64, 10624953561019755799u64, 2789598613442376532u64];
    
    // Create initial state (2 field elements each for x and y)
    // let mut current_x = [F::random(&mut witness_rng), F::random(&mut witness_rng)];
    // let mut current_y = [F::random(&mut witness_rng), F::ZERO]; // Salt is 0
    let mut current_x = [F::from(1u64), F::from(1u64)];
    let mut current_y = [F::from(3u64), F::from(423u64)];



    // fuck, index 0 can't be used for shift
    // TODO: FIX Rotatable FIELD RR
    // 使用 usable_indices 来设置初始状态
    let init_idx = usable_indices[0];
    w1_values[init_idx] = current_x[0];
    w2_values[init_idx] = current_x[1];
    w3_values[init_idx] = current_y[0];
    w4_values[init_idx] = current_y[1];
    wo_values[init_idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    
    // Simple constraint for input gate
    q1[init_idx] = F::ONE;
    q2[init_idx] = F::ONE;
    q3[init_idx] = F::ONE;
    q4[init_idx] = F::ONE;
    qo[init_idx] = F::ONE;

    // Create self-referencing permutations for witness polynomials
    for poly_idx in 14..19 {  // w1, w2, w3, w4, wo
        permutation.copy((poly_idx, init_idx), (poly_idx, init_idx));
    }
    


    // Anemoi rounds - simulate the actual Anemoi permutation
    for round in 0..NUM_ROUNDS {
        let current_gate_idx = usable_indices[round + 1];
        
        // 确保下一个索引也存在（用于 Rotation::next()）
        if round + 2 >= usable_indices.len() {
            panic!("Not enough usable indices for round {}", round);
        }
        
        // Set preprocessed round key selectors (use preprocessed values for qprk)
        qprk1[current_gate_idx] = preprocessed_round_keys_x[round][0];
        qprk2[current_gate_idx] = preprocessed_round_keys_x[round][1];
        qprk3[current_gate_idx] = preprocessed_round_keys_y[round][0];
        qprk4[current_gate_idx] = preprocessed_round_keys_y[round][1];
        
        // Store pre-round values
        w1_values[current_gate_idx] = current_x[0];
        w2_values[current_gate_idx] = current_x[1];
        w3_values[current_gate_idx] = current_y[0];
        w4_values[current_gate_idx] = current_y[1];

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
            current_y[i] -= current_x[i].pow(alpha_inv); // Use alpha_inv for power
            current_x[i] += generator * current_y[i].square() + generator_inv;
        }

        // // 测试anemoi
        // // Verify Anemoi constraints after each round
        // let w1 = w1_values[idx];
        // let w2 = w2_values[idx];
        // let w3 = w3_values[idx];
        // let w4 = w4_values[idx];
        // let qprk1_val = qprk1[idx];
        // let qprk2_val = qprk2[idx];
        // let qprk3_val = qprk3[idx];
        // let qprk4_val = qprk4[idx];
        
        // // Next round values - these should be the state after this round's transformation
        // let w1_next = current_x[0];  // Next round w1
        // let w2_next = current_x[1];  // Next round w2
        // let w3_next = current_y[0];  // Next round w3
        // let wo = current_y[1];

        
        // // wo for the current gate - this should satisfy the output constraint
        // // For now, let's see what value would make the constraints zero

        // // Helper expressions for Anemoi round
        // let c_prime_1 = w1 + w4 + generator * (w2 + w3) + qprk3_val;
        // let c_prime_2 = generator * (w1 + w4) + generator_square_plus_one * (w2 + w3) + qprk4_val;
        
        // // Calculate differences for power operations
        // let diff1 = c_prime_1 - w3_next;
        // let diff2 = c_prime_2 - wo;
        
        // // First Anemoi equation
        // let anemoi_1 = qprk3_val * (
        //     diff1.pow([5u64])  // (c_prime_1 - w3_next)^5
        //     + generator * c_prime_1.square()
        //     - (F::from(2u64) * w1 + w4 + generator * (F::from(2u64) * w2 + w3) + qprk1_val)
        // );
        
        // // Second Anemoi equation  
        // let anemoi_2 = qprk3_val * (
        //     diff2.pow([5u64])  // (c_prime_2 - wo)^5
        //     + generator * c_prime_2.square()
        //     - (generator * (F::from(2u64) * w1 + w4) 
        //         + generator_square_plus_one * (F::from(2u64) * w2 + w3) + qprk2_val)
        // );
        
        // // Third Anemoi equation  
        // let anemoi_3 = qprk3_val * (
        //     diff1.pow([5u64])  // (c_prime_1 - w3_next)^5
        //     + generator * w3_next.square()
        //     + generator_inv
        //     - w1_next
        // );
        
        // // Fourth Anemoi equation
        // let anemoi_4 = qprk3_val * (
        //     diff2.pow([5u64])  // (c_prime_2 - wo)^5
        //     + generator * wo.square()
        //     + generator_inv
        //     - w2_next
        // );

        // // Try to find the correct wo value that would make anemoi_2 and anemoi_4 zero
        // // From anemoi_2: qprk3 * ((c_prime_2 - wo)^5 + g * c_prime_2^2 - (...)) = 0
        // // From anemoi_4: qprk3 * ((c_prime_2 - wo)^5 + g * wo^2 + g_inv - w2_next) = 0
        
        // println!("Round {}: Constraint verification", round);
        // println!("  Pre-round state: w1={:?}, w2={:?}, w3={:?}, w4={:?}", w1, w2, w3, w4);
        // println!("  Post-round state: x={:?}, y={:?}", current_x, current_y);
        // println!("  Next values: w1_next={:?}, w2_next={:?}, w3_next={:?}", w1_next, w2_next, w3_next);
        // println!("  wo (Jive sum)={:?}", wo);
        // println!("  Helper values:");
        // println!("    c_prime_1={:?}", c_prime_1);
        // println!("    c_prime_2={:?}", c_prime_2);
        // println!("    diff1 (c_prime_1 - w3_next)={:?}", diff1);
        // println!("    diff2 (c_prime_2 - wo)={:?}", diff2);
        // println!("  Preprocessed round keys:");
        // println!("    qprk1={:?}", qprk1_val);
        // println!("    qprk2={:?}", qprk2_val);
        // println!("    qprk3={:?}", qprk3_val);
        // println!("    qprk4={:?}", qprk4_val);
        // println!("  CONSTRAINT VALUES:");
        // println!("    anemoi_1 = {:?}", anemoi_1);
        // println!("    anemoi_2 = {:?}", anemoi_2);  
        // println!("    anemoi_3 = {:?}", anemoi_3);
        // println!("    anemoi_4 = {:?}", anemoi_4);
        
        // 5.24 10:03 运算正确
        // println!("current_x: {:?}", current_x);
        // println!("current_y: {:?}", current_y);   
    }

    // Set up next round values for the last round to enable Rotation::next() constraints
    let final_idx = usable_indices[NUM_ROUNDS + 1];
    w1_values[final_idx] = current_x[0];
    w2_values[final_idx] = current_x[1];
    w3_values[final_idx] = current_y[0];
    w4_values[final_idx] = current_y[1];

    // Final output gate
    let output_idx = usable_indices[NUM_ROUNDS + 2];
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

    w1_values[output_idx] = current_x[0];
    w2_values[output_idx] = current_x[1];
    w3_values[output_idx] = current_y[0];
    w4_values[output_idx] = current_y[1];
    wo_values[output_idx] = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    
    // Output constraint
    q1[output_idx] = F::ONE;
    q2[output_idx] = F::ONE;
    q3[output_idx] = F::ONE;
    q4[output_idx] = F::ONE;
    qo[output_idx] = F::ONE;


    // Set up copy constraints using usable indices
    for round in 1..=NUM_ROUNDS {
        let current_idx = usable_indices[round];
        let next_idx = usable_indices[round + 1];
        
        // wo -> w4 connection (Section 5.3 optimization)
        wo_values[current_idx] = w4_values[next_idx];
        permutation.copy((18, current_idx), (17, next_idx)); // wo -> w4 (next round)
    }
    
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
        constraints: anemoi_hash_circuit_info::<F>(num_vars).constraints,
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


/// Generate padding constants for domain separation using π digits
/// Following Section 7.2 of the paper
fn generate_padding_constants<F: PrimeField>(num_levels: usize) -> Vec<F> {
    // In practice, these would be derived from π digits as described in the paper
    // This is a simplified implementation using deterministic constants
    (0..num_levels)
        .map(|i| F::ZERO) // Simple deterministic pattern
        .collect()
}

/// Complete Merkle tree membership proof circuit with full Jive CRH implementation
/// Verifies that a leaf is part of a 3-ary Merkle tree with given root
pub fn merkle_membership_proof_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    leaf_value: F,
    merkle_path: Vec<MerkleProofNode<F>>, // Path from leaf to root
    expected_root: F,
    mut _preprocess_rng: impl RngCore,
    mut _witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    
    const NUM_ROUNDS: usize = 14;
    const GATES_PER_JIVE: usize = NUM_ROUNDS + 4; // input + 14 rounds + sum + output
    const GATES_PER_LEVEL: usize = GATES_PER_JIVE + 10; // Jive CRH + position verification gates
    
    let tree_depth = merkle_path.len();
    let padding_constants = vec![F::ZERO; tree_depth];
    
    // Calculate total gates needed
    let total_gates = tree_depth * GATES_PER_LEVEL + 10; // Buffer for root verification
    let num_vars = (total_gates as f64).log2().ceil() as usize;
    let size = 1 << num_vars;
    let usable_indices = R::from(num_vars).usable_indices();
    
    if usable_indices.len() < total_gates {
        panic!("Not enough usable indices: {} < {}", usable_indices.len(), total_gates);
    }
    
    // Initialize witness polynomials
    let mut w1_values = vec![F::ZERO; size];
    let mut w2_values = vec![F::ZERO; size];
    let mut w3_values = vec![F::ZERO; size];
    let mut w4_values = vec![F::ZERO; size];
    let mut wo_values = vec![F::ZERO; size];
    
    // Initialize selector polynomials
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
    let generator_inv = F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap_or(F::ZERO);
    let generator_square_plus_one = F::from(26u64);
    
    // Real ROUND_KEYS (same as in your paste)
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
    
    let alpha_inv = [14981214993055009997u64, 6006880321387387405u64, 10624953561019755799u64, 2789598613442376532u64];
    
    let mut gate_counter = 0;
    let mut current_hash = leaf_value;
    
    // Process each level of the Merkle path
    for (level, proof_node) in merkle_path.iter().enumerate() {
        let level_padding = padding_constants[level];
        
        // === Position Verification Gates ===
        
        // Gate 1: Boolean constraint for position indicators
        let pos_gate_idx = usable_indices[gate_counter];
        gate_counter += 1;
        
        let (is_left, is_middle, is_right) = match proof_node.position {
            MerklePosition::Left => (F::ONE, F::ZERO, F::ZERO),
            MerklePosition::Middle => (F::ZERO, F::ONE, F::ZERO),
            MerklePosition::Right => (F::ZERO, F::ZERO, F::ONE),
        };
        
        w1_values[pos_gate_idx] = F::ZERO; // Zero variable
        w2_values[pos_gate_idx] = is_left;
        w3_values[pos_gate_idx] = is_middle;
        w4_values[pos_gate_idx] = is_right;
        wo_values[pos_gate_idx] = F::ZERO;
        
        // Constraint: is_left + is_middle + is_right = 1
        q2[pos_gate_idx] = F::ONE;
        q3[pos_gate_idx] = F::ONE;
        q4[pos_gate_idx] = F::ONE;
        qc[pos_gate_idx] = -F::ONE;
        
        // Boolean constraints: each indicator is 0 or 1
        qb[pos_gate_idx] = F::ONE;
        
        // Determine the three children based on position
        let (left_hash, middle_hash, right_hash) = match proof_node.position {
            MerklePosition::Left => (proof_node.current_hash , proof_node.sibling1, proof_node.sibling2),
            MerklePosition::Middle => (proof_node.sibling1, proof_node.current_hash, proof_node.sibling2),
            MerklePosition::Right => (proof_node.sibling1, proof_node.sibling2, proof_node.current_hash),
        };

        // Gate 2: First position selection
        let sel1_gate_idx = usable_indices[gate_counter];
        gate_counter += 1;
        
        w1_values[sel1_gate_idx] = left_hash;
        w2_values[sel1_gate_idx] = is_left;
        w3_values[sel1_gate_idx] = middle_hash;
        w4_values[sel1_gate_idx] = is_middle;
        wo_values[sel1_gate_idx] = left_hash * is_left + middle_hash * is_middle;
        
        qm1[sel1_gate_idx] = F::ONE;
        qm2[sel1_gate_idx] = F::ONE;
        qo[sel1_gate_idx] = F::ONE;
        
        // Gate 3: Second position selection  
        let sel2_gate_idx = usable_indices[gate_counter];
        gate_counter += 1;
        
        let partial_sum = wo_values[sel1_gate_idx];
        w1_values[sel2_gate_idx] = right_hash;
        w2_values[sel2_gate_idx] = is_right;
        w3_values[sel2_gate_idx] = partial_sum;
        w4_values[sel2_gate_idx] = F::ZERO;
        wo_values[sel2_gate_idx] = current_hash;

        // println!("current_hash: {:?}", current_hash);
        // println!("partial_sum: {:?}", partial_sum);
        // println!("right_hash: {:?}", right_hash);
        
        qm1[sel2_gate_idx] = F::ONE;
        q3[sel2_gate_idx] = F::ONE;
        qo[sel2_gate_idx] = F::ONE;
        
        // === Complete Jive CRH Implementation ===
        
        // Input gate for Jive CRH
        let jive_input_idx = usable_indices[gate_counter];
        gate_counter += 1;
        // let jive_start = usable_indices[gate_counter];
        
        let mut jive_x = [left_hash, middle_hash];
        let mut jive_y = [right_hash, level_padding];
        let sum_before_perm = jive_x[0] + jive_x[1] + jive_y[0] + jive_y[1];
        
        w1_values[jive_input_idx] = jive_x[0];
        w2_values[jive_input_idx] = jive_x[1];
        w3_values[jive_input_idx] = jive_y[0];
        w4_values[jive_input_idx] = jive_y[1];
        wo_values[jive_input_idx] = sum_before_perm;
        
        // Enforce padding constant constraint
        q1[jive_input_idx] = F::ONE;
        q2[jive_input_idx] = F::ONE;
        q3[jive_input_idx] = F::ONE;
        q4[jive_input_idx] = F::ONE;
        qo[jive_input_idx] = F::ONE;
        
        // q4[jive_input_idx] = F::ONE;
        // qc[jive_input_idx] = -level_padding;
        
        // Create self-referencing permutations for this Jive instance
        for poly_idx in 14..19 {
            permutation.copy((poly_idx, jive_input_idx), (poly_idx, jive_input_idx));
        }
        
        // 14 Anemoi rounds
        for round in 0..NUM_ROUNDS {
            let round_gate_idx = usable_indices[gate_counter];
            gate_counter += 1;
            
            // Set preprocessed round keys
            qprk1[round_gate_idx] = preprocessed_round_keys_x[round][0];
            qprk2[round_gate_idx] = preprocessed_round_keys_x[round][1];
            qprk3[round_gate_idx] = preprocessed_round_keys_y[round][0];
            qprk4[round_gate_idx] = preprocessed_round_keys_y[round][1];
            
            // Store pre-round state
            w1_values[round_gate_idx] = jive_x[0];
            w2_values[round_gate_idx] = jive_x[1];
            w3_values[round_gate_idx] = jive_y[0];
            w4_values[round_gate_idx] = jive_y[1];
            
            // Apply Anemoi round transformation
            jive_x[0] += round_keys_x[round][0];
            jive_x[1] += round_keys_x[round][1];
            jive_y[0] += round_keys_y[round][0];
            jive_y[1] += round_keys_y[round][1];
            
            // MDS matrix
            let temp_x0 = jive_x[0] + generator * jive_x[1];
            let temp_x1 = generator * jive_x[0] + generator_square_plus_one * jive_x[1];
            jive_x[0] = temp_x0;
            jive_x[1] = temp_x1;
            
            let temp_y0 = jive_y[1] + generator * jive_y[0];
            let temp_y1 = generator * jive_y[1] + generator_square_plus_one * jive_y[0];
            jive_y[0] = temp_y0;
            jive_y[1] = temp_y1;
            
            // PHT transformation
            for i in 0..2 {
                jive_y[i] += jive_x[i];
                jive_x[i] += jive_y[i];
            }
            
            // S-box
            for i in 0..2 {
                jive_x[i] -= generator * jive_y[i].square();
                jive_y[i] -= jive_x[i].pow(alpha_inv);
                jive_x[i] += generator * jive_y[i].square() + generator_inv;
            }
        }
        
        // Next round values gate for rotation::next() constraints
        let next_round_idx = usable_indices[gate_counter];
        gate_counter += 1;
        
        w1_values[next_round_idx] = jive_x[0];
        w2_values[next_round_idx] = jive_x[1];
        w3_values[next_round_idx] = jive_y[0];
        w4_values[next_round_idx] = jive_y[1];

        // final state constraint
        // TODO: 这里为了安全性最好做一下拆分，目前是PoC，所以直接写在一起了
        q1[next_round_idx] = F::from(3u64) * generator + F::from(3u64);
        q2[next_round_idx] = F::from(3u64) * (generator * generator + generator + F::ONE);
        q3[next_round_idx] = F::from(2u64) * (generator * generator + generator + F::ONE);
        q4[next_round_idx] = F::from(2u64) * generator + F::from(2u64);
        qo[next_round_idx] = F::ONE;
  
        // Final MDS and PHT
        let temp_x0 = jive_x[0] + generator * jive_x[1];
        let temp_x1 = generator * jive_x[0] + generator_square_plus_one * jive_x[1];
        jive_x[0] = temp_x0;
        jive_x[1] = temp_x1;
        
        let temp_y0 = jive_y[1] + generator * jive_y[0];
        let temp_y1 = generator * jive_y[1] + generator_square_plus_one * jive_y[0];
        jive_y[0] = temp_y0;
        jive_y[1] = temp_y1;
        
        for i in 0..2 {
            jive_y[i] += jive_x[i];
            jive_x[i] += jive_y[i];
        }
        
        let sum_after_perm = jive_x[0] + jive_x[1] + jive_y[0] + jive_y[1];
        wo_values[next_round_idx] = sum_after_perm;
        
        
        // Final Jive output gate
        let jive_output_idx = usable_indices[gate_counter];
        gate_counter += 1;
        
        current_hash = sum_before_perm + sum_after_perm; // Jive output
        
        w1_values[jive_output_idx] = sum_before_perm;
        w2_values[jive_output_idx] = sum_after_perm;
        w3_values[jive_output_idx] = F::ZERO;
        w4_values[jive_output_idx] = F::ZERO;
        wo_values[jive_output_idx] = current_hash;
        // println!("current_hash: {:?}", current_hash);
        
        // Constraint: wo = w1 + w2
        q1[jive_output_idx] = F::ONE;
        q2[jive_output_idx] = F::ONE;
        qo[jive_output_idx] = F::ONE;

        permutation.copy((18, next_round_idx), (15, jive_output_idx));
        permutation.copy((14, jive_output_idx), (18, usable_indices[gate_counter - 17]));
        
        // Set up copy constraints for this Jive instance  
        // for round in 1..=NUM_ROUNDS {
        //     let current_idx = usable_indices[jive_input_idx + round];
        //     let next_idx = usable_indices[jive_input_idx + round + 1];
            
        //     wo_values[current_idx] = w4_values[next_idx];
        //     permutation.copy((18, current_idx), (17, next_idx)); // wo -> w4
        // }
    }
    
    // Final root verification gate
    let root_gate_idx = usable_indices[gate_counter];
    w1_values[root_gate_idx] = current_hash;
    w2_values[root_gate_idx] = expected_root;
    w3_values[root_gate_idx] = F::ZERO;
    w4_values[root_gate_idx] = F::ZERO;
    wo_values[root_gate_idx] = F::ZERO;
    
    // Constraint: current_hash - expected_root = 0
    q1[root_gate_idx] = F::ONE;
    q2[root_gate_idx] = -F::ONE;


    
    let circuit_info = PlonkishCircuitInfo {
        k: num_vars,
        num_instances: vec![],
        preprocess_polys: vec![
            q1, q2, q3, q4, qo, qm1, qm2, qc, qecc, qb,
            qprk1, qprk2, qprk3, qprk4
        ],
        num_witness_polys: vec![5],
        num_challenges: vec![0],
        constraints: anemoi_hash_circuit_info_original::<F>(num_vars).constraints,
        lookups: Vec::new(),
        permutations: permutation.into_cycles(),
        max_degree: Some(7),
    };
    
    let witness = vec![w1_values, w2_values, w3_values, w4_values, wo_values];
    
    (circuit_info, MockCircuit::new(vec![], witness))
}




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

/// Generate random Merkle membership proof test data using generic PrimeField
/// Returns (leaf_value, merkle_path, expected_root)
pub fn generate_random_merkle_proof<F: PrimeField>(
    depth: usize,
    mut rng: impl RngCore,
) -> (F, Vec<MerkleProofNode<F>>, F) {
    // assert!(depth > 0, "Merkle tree depth must be positive");
    // assert!(depth <= 32, "Maximum practical depth is 32");
    
    // Generate padding constants for domain separation
    let padding_constants = vec![F::ZERO; depth];
    
    // Generate random leaf value
    let leaf_value = F::random(&mut rng);
    
    // Generate random path (position at each level from leaf to root)
    let path_positions: Vec<MerklePosition> = (0..depth)
        .map(|_| match rng.next_u32() % 3 {
            0 => MerklePosition::Left,
            1 => MerklePosition::Middle,
            _ => MerklePosition::Right,
        })
        .collect();
    
    // Build merkle path from leaf to root
    let mut current_hash = leaf_value;
    let mut merkle_path = Vec::new();
    
    for (level, position) in path_positions.iter().enumerate() {
        // Generate random siblings
        let sibling1 = F::random(&mut rng);
        let sibling2 = F::random(&mut rng);
        
        // Determine the three children based on current node's position
        let (left_hash, middle_hash, right_hash) = match position {
            MerklePosition::Left => (current_hash, sibling1, sibling2),
            MerklePosition::Middle => (sibling1, current_hash, sibling2),
            MerklePosition::Right => (sibling1, sibling2, current_hash),
        };
        
        // Use level-specific padding constant for domain separation
        let level_padding = padding_constants[level];
        
        // Compute parent hash using Anemoi Jive CRH
        // This follows the 3-ary Merkle tree construction in Section 7.1
        let parent_hash = {
            let x = [left_hash, middle_hash];  // First 2 children
            let y = [right_hash, level_padding]; // Third child + domain separation
            anemoi_jive_hash(x, y)
        };
        
        // Add to merkle path (this represents the current level)
        merkle_path.push(MerkleProofNode {
            position: position.clone(),
            current_hash,  // Hash of current node (before going up to parent)
            sibling1,      // First sibling
            sibling2,      // Second sibling
        });
        
        // Move up to parent for next iteration
        current_hash = parent_hash;
    }
    
    let expected_root = current_hash;
    
    // hash 计算没有问题
    // println!("leaf_value: {:?}", leaf_value);
    // println!("merkle_path: {:?}", merkle_path);
    // println!("expected_root: {:?}", expected_root);

    (leaf_value, merkle_path, expected_root)
}



/// Generate random Merkle membership proof circuit with proper generic implementation
pub fn rand_merkle_membership_proof_circuit<F: PrimeField, R: Rotatable + From<usize>>(
    depth: usize,
    mut preprocess_rng: impl RngCore,  
    mut witness_rng: impl RngCore,
) -> (PlonkishCircuitInfo<F>, impl PlonkishCircuit<F>) {
    // Generate test data using generic implementation
    let (leaf_value, merkle_path, expected_root) = 
        generate_random_merkle_proof(depth, &mut witness_rng);
    
    merkle_membership_proof_circuit::<F, R>(
        leaf_value,
        merkle_path,
        expected_root,
        preprocess_rng,
        witness_rng,
    )
}


/// Generic Anemoi Jive hash implementation using PrimeField
/// This extracts the logic from merkle_membership_proof_circuit to make it reusable
pub fn anemoi_jive_hash<F: PrimeField>(x: [F; 2], y: [F; 2]) -> F {
    const NUM_ROUNDS: usize = 14;
    
    // Use real Anemoi constants (converted to generic field)
    let generator = F::from(5u64);
    let generator_inv = F::from_str_vartime("8755297148735710088898562298102910035419345760166413737479281674630323398247").unwrap_or(F::ZERO);
    let generator_square_plus_one = F::from(26u64);
    
    // Real ROUND_KEYS_X and ROUND_KEYS_Y from AnemoiJive256
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
    
    // Alpha inverse values for S-box
    let alpha_inv = [14981214993055009997u64, 6006880321387387405u64, 10624953561019755799u64, 2789598613442376532u64];
    
    // Set initial state
    let mut current_x = x;
    let mut current_y = y;
    let sum_before_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    
    // Anemoi rounds - simulate the actual Anemoi permutation
    for round in 0..NUM_ROUNDS {
        // Apply Anemoi round transformation following the implementation
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
        
        // Step 3: PHT transformation
        for i in 0..2 {
            current_y[i] += current_x[i];
            current_x[i] += current_y[i];
        }
        
        // Step 4: S-box application
        for i in 0..2 {
            current_x[i] -= generator * current_y[i].square();
            current_y[i] -= current_x[i].pow(alpha_inv); // Use alpha_inv for power
            current_x[i] += generator * current_y[i].square() + generator_inv;
        }
    }
    
    // Final transformations (post-rounds)
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
    
    let sum_after_perm = current_x[0] + current_x[1] + current_y[0] + current_y[1];
    
    // Jive output: sum of input and output of permutation
    sum_before_perm + sum_after_perm
}



/// Helper structures for Merkle tree proofs
#[derive(Debug, Clone)]
pub struct MerkleProofNode<F: PrimeField> {
    pub position: MerklePosition,
    pub current_hash: F, // Current node hash
    pub sibling1: F, // Left or first sibling
    pub sibling2: F, // Right or second sibling  
}

#[derive(Debug, Clone)]
pub enum MerklePosition {
    Left,   // Current node is left child
    Middle, // Current node is middle child
    Right,  // Current node is right child
}
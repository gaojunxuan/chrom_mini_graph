use fxhash::FxHashMap;
use mini_alph::chain;
use mini_alph::seeding_methods_bit;
use mini_alph::simulation_utils_bit;

fn main() {
    let theta = 0.05;
    let n = 1000;
    let w = 9;
    let k = 16;
    let s = simulation_utils_bit::gen_rand_string(n);
    let printable: Vec<u8> = s.iter().collect();
    let (s_p, mutated_pos) = simulation_utils_bit::gen_mutated_string(&s, theta);
    let printable_p: Vec<u8> = s_p.iter().collect();

    //    println!("{}",s.to_string());
    //    println!("{:?}",printable_p);

//        let (seeds1, positions_selected) = seeding_methods_bit::minimizer_seeds(&s, w, k);
    let (seeds1, positions_selected) = seeding_methods_bit::open_sync_seeds(&s, k, 3);
//        let (seeds2, positions_selected_p) = seeding_methods_bit::minimizer_seeds(&s_p, w, k);
    let (seeds2, positions_selected_p) = seeding_methods_bit::open_sync_seeds(&s_p, k, 3);

    //    dbg!(&positions_selected);

    chain::chain_seeds(
        &seeds1,
        &seeds2,
        &positions_selected,
        &positions_selected_p,
        w,
        30,
        k,
    );

    let mut list_mutpos = vec![];
    for i in 0..mutated_pos.len() {
        if !mutated_pos[i] {
            list_mutpos.push(i)
        }
    }

    println!("{:?}", list_mutpos);

    //    println!("{:?}",minimizers);
    //    println!("{:?}",positions_selected);
    //    dbg!(minimizers.len());
    //    let mut distance_vec = vec![0; w + 1];
    //
    //    for i in 0..positions_selected.len() - 1 {
    //        distance_vec[positions_selected[i + 1] - positions_selected[i]] += 1;
    //    }
    //
    //    dbg!(distance_vec);
}

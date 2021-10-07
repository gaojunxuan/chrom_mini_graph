use rand::{thread_rng};
use rand::distributions::{Bernoulli, Distribution, Uniform};
use debruijn::dna_string::*;


pub fn gen_rand_string(n: usize) ->  DnaString{

    let mut rng = rand::thread_rng();
    let rv = Uniform::from(0..4);
    
    let mut return_string = DnaString::new();
    for _i in 0..n{
        return_string.push(rv.sample(&mut rng));
    }
    return return_string;
}

pub fn gen_mutated_string(sequence: &DnaString, theta: f64) -> (DnaString, Vec<bool>) {
    let mut return_vec = DnaString::new();
    let mut return_vec_bool = vec![];
    let d = Bernoulli::new(theta).unwrap();
    let mut rng = thread_rng();
    //Turn A to T, and C to G with probability theta.
    for base in sequence.iter() {
        let x = d.sample(&mut rng);
        let mut new_base = 5;
        if x {
            if base == 0 {
                new_base = 3;
            }
            if base == 1 {
                new_base = 2;
            }
            if base == 3 {
                new_base = 0;
            }
            if base == 2 {
                new_base = 1;
            }

            return_vec.push(new_base);
            return_vec_bool.push(false);
        } else {
            return_vec.push(base);
            return_vec_bool.push(true);
        }
    }
    return (return_vec, return_vec_bool);
}


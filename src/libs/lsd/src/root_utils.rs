use cmg_shared::data_structs::{ColoredGraph, ColoredCmg};
use crate::{constants::WHITE, cycle_utils::{find_cover_cut, CycleGenerator}, path_utils::find_paths};

pub struct RootGenerator {
    next_id: Option<u32>,
    max_id: u32,
    cycle_generator: CycleGenerator,
}

impl RootGenerator {
    pub fn new(n: u32) -> RootGenerator {
        match n {
            0 => {
                RootGenerator {
                    next_id: None,
                    max_id: 0,
                    cycle_generator: CycleGenerator::new(0),
                }
            },
            _ => {
                RootGenerator {
                    next_id: Some(0),
                    max_id: (n - 1) as u32,
                    cycle_generator: CycleGenerator::new(n),
                }
            }
        }
    }

    pub fn generate_next_root(&mut self, g: &mut ColoredCmg) -> Option<(u32, bool)> {
        // println!("Generating next root. Next ID: {:?}, Max ID: {:?}", self.next_id, self.max_id);
        while self.next_id.is_some() {
            if let Some(v) = self.next_id {
                if v < self.max_id {
                    self.next_id = Some(v + 1);
                } else {
                    self.next_id = None;
                }
                // println!("Current node color is {:?}, in_degree is {:?}", g.get_color(v), g.graph.in_degree(v));
                if g.get_color(v) == WHITE {
                    if g.graph.in_degree(v) == 0 {
                        return Some((v as u32, false));
                    }
                }
            }
        }
        
        while let Some(cycle) = self.cycle_generator.find_next_cycle(g) {
            let converted_cycle = cycle.iter().map(|x| *x).collect::<Vec<_>>();
            let (cover, value) = find_paths(g, &converted_cycle);
            if cover {
                return Some((cycle[*value.first().unwrap() as usize], true));
            } else {
                let paths = value;
                let (cutpoint, position) = find_cover_cut(&paths);
                if cutpoint {
                    return Some((cycle[position as usize], true));
                } else {
                    return Some((cycle[position as usize], false));
                }
            }
        }
        None
    }
}

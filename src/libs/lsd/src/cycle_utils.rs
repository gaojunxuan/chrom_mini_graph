use std::future::Future;
use genawaiter::{sync::{gen, Gen}, yield_};
use cmg_shared::data_structs::{ColoredGraph, ColoredCmg};
use crate::{constants::{WHITE, GREEN, YELLOW}, hashed_list::HashedList};

/// Find cycles in a graph
/// 
/// `CycleGenerator` itself does not borrow the graph, but the `find_next_cycle` method does.
/// 
pub struct CycleGenerator {
    next_id: Option<u32>,
    max_id: u32,
    stack: Vec<u32>,
}

impl CycleGenerator {
    /// Create a new instance of `CycleGenerator`
    /// 
    /// # Arguments
    /// - `n`: The number of nodes in the graph
    pub fn new(n: u32) -> CycleGenerator {
        match n {
            0 => {
                CycleGenerator {
                    next_id: None,
                    max_id: 0,
                    stack: vec![],
                }
            },
            _ => {
                CycleGenerator {
                    next_id: Some(0),
                    max_id: (n - 1) as u32,
                    stack: vec![],
                }
            }
        }
    }

    pub fn find_next_cycle(&mut self, g: &mut ColoredCmg) -> Option<HashedList<u32>> {
        while self.next_id.is_some() {
            if let Some(node_id) = self.next_id {
                if g.get_color(node_id) == WHITE {
                    while let Some(u) = self.stack.pop() {
                        let c = g.get_color(u);
                        if c == WHITE {
                            g.set_color(u, GREEN);
                            for w in g.graph.get_successors(u) {
                                let wc = g.get_color(w);
                                if wc == WHITE {
                                    self.stack.push(w);
                                } else if wc == GREEN {
                                    let mut cycle = vec![];
                                    while let Some(x) = self.stack.pop() {
                                        if x != w {
                                            if g.get_color(x) == GREEN {
                                                cycle.push(x);
                                            }
                                        }
                                    }
                                    cycle.push(self.stack.pop().unwrap());
                                    // yield reversed(cycle)
                                    return Some(HashedList::new(Some(cycle.iter().rev().map(|x| *x as u32).collect::<Vec<u32>>())));
                                }
                            }
                        } else {
                            if c == GREEN {
                                g.set_color(u, YELLOW);
                            }
                            self.stack.pop();
                        }
                    }
                }
                if node_id < self.max_id {
                    self.next_id = Some(node_id + 1);
                } else {
                    self.next_id = None;
                }
            }
        }
        return None;
    }
}

// // LSD/cycles
// pub fn find_cycles<'a>(g: &'a mut ColoredGraph) -> Gen<Vec<u32>, (), impl Future<Output = ()> + 'a> {
//     return gen!({
//         for node_id in 0..g.graph.nodes.len() {
            
//         }
//     });
// }

pub fn cycle_distance(k: u32, i: u32, end: u32) -> u32 {
    if end > i {
        return end - i - 1;
    }
    return end + k - i - 1;
}

pub fn generate_cycle_range(k: u32, pos: u32, end: u32) -> Gen<u32, (), impl Future<Output = ()>> {
    return gen!({
        let mut i = pos;
        if end < pos {
            while i < end || i >= pos {
                yield_!(i);
                i = (i + 1) % k;
            }
        }
        else {
            while i < end {
                yield_!(i);
                i += 1;
            }
        }
    });
}

// LSD/cycles/cover.py
pub fn find_cover_cut(cmax: &[u32]) -> (bool, u32) {
    let k = cmax.len() as u32;
    
    // Calculate maximal global C-Path
    let m = (0..k)
        .max_by_key(|&i| cycle_distance(k, i as u32, cmax[i as usize]))
        .unwrap();
    let md = cycle_distance(k, m, cmax[m as usize]);
    let mut laststart = m;
    let mut lastend = m;
    
    // While laststart end not in maximal global C-Path
    while md <= cycle_distance(k, m, cmax[laststart as usize]) {
        // Calc maximal overtower
        let (mut next_path, mut max_dist) = (None, 0);
        
        for i in generate_cycle_range(k, lastend, cmax[laststart as usize]) {
            // Filter included paths
            if cycle_distance(k, i, cmax[i as usize]) <= cycle_distance(k, i, cmax[laststart as usize]) {
                continue;
            }
            
            // Calc overtower and save max
            let dist = cycle_distance(k, cmax[laststart as usize], cmax[i as usize]);
            if dist > max_dist {
                max_dist = dist;
                next_path = Some(i);
            }
        }
        
        // If nothing overtower
        if let Some(path) = next_path {
            // Set next path to check
            lastend = cmax[laststart as usize];
            laststart = path;
        } else {
            // Return cut point
            return (true, cmax[laststart as usize]);
        }
    }
    
    // Return Cycle cover
    (false, cmax[laststart as usize])
}
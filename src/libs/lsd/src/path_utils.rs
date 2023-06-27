use std::vec;
use cmg_shared::data_structs::{ColoredGraph, ColoredCmg};
use crate::{constants::{RED, ORANGE, BLACK}, cycle_utils::cycle_distance, dfs_tree::create_tree, hashed_list::HashedList};

fn finish_subtree(g: &mut ColoredCmg, v: u32) {
    let minval = g.graph.get_property(v, "min".to_string()).unwrap_or(u32::MIN);
    let maxval = g.graph.get_property(v, "max".to_string()).unwrap_or(u32::MAX);
    let mut nexts = vec![v];
    g.set_color(v, RED);
    
    while let Some(u) = nexts.pop() {
        g.graph.set_property(u, "min".to_string(), minval);
        g.graph.set_property(u, "max".to_string(), maxval);
        
        for w in g.graph.get_successors(u).clone() {
            if g.get_color(w) == ORANGE {
                g.set_color(w, RED);
                nexts.push(w);
            }
        }
    }
}

fn calculate_cmax(g: &mut ColoredCmg, cycle: &[u32], postorder: &HashedList<u32>, inorder: &HashedList<u32>, i: u32) -> (bool, Option<u32>) {
    let k = cycle.len() as u32;

    let cycle_min = |args: Vec<u32>| -> Option<u32> {
        let mut min_value = u32::MAX;
        let mut min_position = None;
        
        for x in args {
            let new_value = cycle_distance(k, i as u32, x);
            if new_value < min_value {
                min_value = new_value;
                min_position = Some(x);
            }
        }
        
        min_position
    };

    let cycle_max = |args: Vec<u32>| -> Option<u32> {
        let mut max_value: i32 = i32::MIN;
        let mut max_position = None;
        
        for x in args {
            let new_value = cycle_distance(k, i as u32, x) as i32;
            if new_value > max_value {
                max_value = new_value;
                max_position = Some(x);
            }
        }
        
        max_position
    };

    let check_backedge = |w: u32, u: u32, postorder: &HashedList<u32>| -> bool {
        postorder.iter().position(|&r| r == w) < postorder.iter().position(|&r| r == u)
    };

    for j in 0..postorder.len() {
        let v = postorder[j];
        let inpos = inorder.iter().position(|&r| r == v);
        g.graph.set_property(v, "link".to_string(), inpos.unwrap() as u32);

        for child in g.graph.get_successors(v).clone() {
            let color = g.get_color(child);
            
            // Is in cycle
            if cycle.contains(&child) {
                let pos = cycle.iter().position(|&x| x == child).unwrap();
                g.graph.update_property(v, "min".to_string(), cycle_min, vec![pos as u32]);
                g.graph.update_property(v, "max".to_string(), cycle_max, vec![pos as u32]);
            }
            // Do nothing when in other run complete finished
            else if color != BLACK {
                // Is back edge
                if color != RED && check_backedge(v, child, &postorder) {
                    let index = inorder.iter().position(|&r| r == child);
                    match index {
                        Some(index) => {
                            g.graph.update_property(v, "link".to_string(), |v: Vec<u32>| v.iter().min().cloned(), vec![index as u32]);
                        },
                        None => {
                            g.graph.update_property(v, "link".to_string(), |v: Vec<u32>| v.iter().min().cloned(), vec![]);
                        }
                    }
                }
                // Is no back edge
                else {
                    // Look for one vertex cover
                    if let Some(child_min) = g.graph.get_property(child, "min".to_string()) {
                        let child_max = g.graph.get_property(child, "max".to_string()).unwrap();
                        let dist_min = cycle_distance(k, i as u32, child_min);
                        let dist_max = cycle_distance(k, i as u32, child_max);
                        if dist_min > dist_max {
                            return (true, Some(child_min));
                        }
                    }

                    // Update values
                    g.graph.update_property(v, "min".to_string(), cycle_min, vec![g.graph.get_property(child, "min".to_string()).unwrap()]);
                    g.graph.update_property(v, "max".to_string(), cycle_max, vec![g.graph.get_property(child, "max".to_string()).unwrap()]);

                    // If not finished
                    if color != RED {
                        g.graph.update_property(v, "link".to_string(), |v: Vec<u32>| v.iter().min().cloned(), vec![g.graph.get_property(child, "link".to_string()).unwrap()]);
                    }
                }
            }
        }

        // Finished subtree
        if g.graph.get_property(v, "link".to_string()).is_some() && g.graph.get_property(v, "link".to_string()).unwrap() == inpos.unwrap() as u32 {
            finish_subtree(g, v);
        }
    }

    (false, g.graph.get_property(cycle[i as usize], "max".to_string()))
}

pub fn find_paths(g: &mut ColoredCmg, cycle: &[u32]) -> (bool, Vec<u32>) {
    let mut result = Vec::new();
    for i in 0..cycle.len() {
        let c = cycle[i];
        let (postorder, inorder) = create_tree(g, c, cycle.to_vec());
        let (cover, value) = calculate_cmax(g, cycle, &postorder, &inorder, i as u32);
        if cover {
            return (true, vec![value.unwrap()]);
        }
        result.push(value.unwrap());
    }
    (false, result)
}
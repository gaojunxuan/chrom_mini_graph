use cmg_shared::data_structs::{ColoredGraph, Graph, GraphNode, ColoredCmg};
use crate::{constants::{ORANGE, BLACK, RED, GRAY}, hashed_list::HashedList };

pub fn create_tree(g: &mut ColoredCmg, c: u32, cycle: Vec<u32>) -> (HashedList<u32>, HashedList<u32>) {
    // let mut postorder = Vec::new();
    // let mut inorder = Vec::new();
    let mut postorder = HashedList::<u32>::new(Some(vec![]));
    let mut inorder = HashedList::<u32>::new(Some(vec![]));
    
    // Add cycle vertex as first in preorder
    inorder.append(c);
    let mut stack = vec![c];
    
    // Create tree
    while let Some(&top) = stack.last() {
        let color = g.get_color(top);
        
        if color != ORANGE {
            g.set_color(top, ORANGE);
            
            // New vertex found
            inorder.append(top);
            
            for child in g.graph.get_successors(top) {
                // No backedge or finished vertex or cycle vertex
                if !matches!(g.get_color(child), BLACK | ORANGE | RED) && !cycle.contains(&child) {
                    stack.push(child);
                }
            }
        } else {
            // Is already orange and so finished now
            // When not already finished
            if !postorder.contains(&top) {
                // Add postorder
                postorder.append(top);
            }
            // Remove from stack
            stack.pop();
        }
    }
    (postorder, inorder)
}

// LSD/dfs_tree/order.py
pub fn create_dfs_order_cycle(v: u32, g: &mut ColoredCmg) -> (HashedList<u32>, Option<u32>) {
    let mut order = HashedList::<String>::new(Some(vec![]));
    // let mut order = Vec::new();
    let mut stack = vec![v];
    let mut cycle = true;
    let mut v2 = None;
    
    while let Some(&u) = stack.last() {
        let c = g.get_color(u);
        
        if c != GRAY && c != BLACK {
            g.set_color(u, GRAY);
            
            for w in g.graph.get_successors(u) {
                if cycle && w == v {
                    v2 = Some(format!("{}_split_duplicate", v));
                    order.append(v2.as_ref().unwrap().clone());
                    cycle = false;
                }
                
                let color = g.get_color(w);
                
                if color != GRAY && color != BLACK {
                    stack.push(w as u32);
                }
            }
        } else {
            stack.pop();
            
            if c == GRAY {
                order.append(u.to_string());
                g.set_color(u, BLACK);
            }
        }
    }
    // remove split duplicates (elements that are not convertible to u32) and convert to u32
    let order = order.iter().filter_map(|x| x.parse::<u32>().ok()).collect::<HashedList<u32>>();
    // convert v2 to u32 after removing _split_duplicate tag
    let v2 = v2.map(|x| x.replace("_split_duplicate", "").parse::<u32>().unwrap());
    (order, v2)
}

pub fn create_dfs_order(v: u32, g: &mut ColoredCmg) -> HashedList<u32> {
    let mut order = HashedList::<u32>::new(Some(vec![]));
    let mut stack = vec![v];
    
    while let Some(&u) = stack.last() {
        let c = g.get_color(u);
        
        if c != GRAY && c != BLACK {
            g.set_color(u, GRAY);
            
            for w in g.graph.get_successors(u) {
                let color = g.get_color(w);
                
                if color != GRAY && color != BLACK {
                    stack.push(w as u32);
                }
            }
        } else {
            stack.pop();
            
            if c == GRAY {
                order.append(u);
                g.set_color(u, BLACK);
            }
        }
    }
    order
}
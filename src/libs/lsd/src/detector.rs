use cmg_shared::data_structs::{ColoredGraph, Graph, Bubble, Cmg, ColoredCmg};
use crate::{dfs_tree::{create_dfs_order_cycle, create_dfs_order}, root_utils::RootGenerator, hashed_list::HashedList};

pub fn detect<'a>(graph: &'a mut Cmg<'a>, bubbles: &mut Vec<Bubble>) {
    let now = std::time::Instant::now();
    let mut g = ColoredCmg::new::<'a>(graph);
    let mut root_generator = RootGenerator::new(g.graph.nodes.len() as u32);
    // let mut bubble_counter = 0;
    let reporter = |x: &[u32], bubbles: &mut Vec<Bubble>| {
        if x.len() > 10 {
            let bubble = Bubble {
                id: bubbles.len(),
                // kmers: x.iter().map(|&x| x as u32).collect::<Vec<u32>>()
                kmers: vec![0],
                start: x.first().unwrap().clone() as u32,
                end: x.last().unwrap().clone() as u32,
                longest_path_length: None,
                shortest_path_length: None,
            };
            bubbles.push(bubble);
        }
    };
    while let Some((v, cycle)) = root_generator.generate_next_root(&mut g) {
        if cycle {
            let (order, v2) = create_dfs_order_cycle(v, &mut g);
            // let order = order.iter().map(|x| *x as usize).collect::<HashedList<_>>();
            // wrap reporter in SCC filter
            superbubble(&g, &order, reporter, v, v2, bubbles);
        } else {
            let order = create_dfs_order(v, &mut g);
            // let order = order.iter().map(|x| *x as usize).collect::<HashedList<_>>();
            superbubble(&g, &order, reporter, v, None, bubbles);
        }
    }
    println!("Found {} bubbles in {} seconds", bubbles.len(), now.elapsed().as_secs_f64());
}

fn superbubble<T>(
    g: &ColoredCmg,
    order: &HashedList<u32>,
    mut reporter: T,
    v: u32,
    v2: Option<u32>,
    bubbles: &mut Vec<Bubble>
)
where
    T: FnMut(&[u32], &mut Vec<Bubble>),
{
    let mut stack: Vec<Option<usize>> = Vec::new();
    let mut out_parent_map: Vec<isize> = Vec::new();
    let mut t: Option<usize> = None;

    for k in 0..order.len() {
        let child = out_child(k, g, order, v, v2);
        
        if child == (k as isize - 1) as isize {
            stack.push(t);
            t = Some(k - 1);
        } else {
            while let Some(t_value) = t {
                if t_value as isize > child {
                    let t2 = stack.pop().unwrap_or(None);
                    if let Some(t2_value) = t2 {
                        out_parent_map[t2_value] = out_parent_map[t_value].max(out_parent_map[t2_value]);
                    }
                    t = t2;
                } else {
                    break;
                }
            }
        }
        
        if let Some(t_value) = t {
            if out_parent_map[t_value] == k as isize {
                let slice = &order[(t_value..=k)].iter().rev().cloned().collect::<Vec<u32>>();
                reporter(slice, bubbles);
                let t2 = stack.pop().unwrap_or(None);
                if let Some(t2_value) = t2 {
                    out_parent_map[t2_value] = out_parent_map[t_value].max(out_parent_map[t2_value]);
                }
                t = t2;
            }
        }
        
        out_parent_map.push(out_parent(k, g, order, v, v2));
        
        if let Some(t_value) = t {
            out_parent_map[t_value] = out_parent_map[t_value].max(out_parent_map[k]);
        }
    }
}

fn out_parent(k: usize, g: &ColoredCmg, order: &HashedList<u32>, v: u32, v2: Option<u32>) -> isize {
    let u = order[k];
    let u = if let Some(v2_value) = v2 {
        if u == v2_value {
            v
        } else {
            u
        }
    } else {
        u
    };
    
    if g.graph.in_degree(u) == 0 {
        return isize::MAX;
    }
    
    let mut maximum = 0;
    
    for w in g.graph.get_predecessors(u) {
        // let pos = order.iter().position(|&x| x == *w).unwrap();
        let pos = order.index_of(&w);
        
        if pos <= k as isize {
            return isize::MAX;
        }
        
        maximum = maximum.max(pos);
    }
    
    maximum
}

fn out_child(k: usize, g: &ColoredCmg, order: &HashedList<u32>, v: u32, v2: Option<u32>) -> isize {
    let u = order[k];
    let u = if let Some(v2_value) = v2 {
        if u == v2_value {
            v
        } else {
            u
        }
    } else {
        u
    };
    
    if g.graph.out_degree(u) == 0 {
        return -2;
    }
    
    let mut minimum = isize::MAX;
    
    for w in g.graph.get_successors(u) {
        let mut pos: isize = 0;
        if let Some(v2_value) = v2 {
            if w == v && v2.is_some() {
                // order.iter().position(|&x| x == v2_value as usize).unwrap()
                pos = order.index_of(&(v2_value));
            } else {
                // order.iter().position(|&x| x == *w).unwrap()
                pos = order.index_of(&w);
            }
        } else {
            // order.iter().position(|&x| x == *w).unwrap()
            pos = order.index_of(&w);
        };
        
        if pos >= k as isize {
            return -2;
        }
        
        minimum = minimum.min(pos as isize);
    }
    
    minimum
}
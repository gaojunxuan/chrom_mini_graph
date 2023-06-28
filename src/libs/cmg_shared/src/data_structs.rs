use std::collections::{HashSet, BTreeSet, BTreeMap, HashMap};
use std::vec::Vec;
use debruijn::kmer::Kmer16;
use smallvec::SmallVec;
use serde::{Serialize, Deserialize};

use crate::constants::WHITE;

//First is ref, second is query
pub type Anchors = Vec<(u32, u32)>;
pub type Color = u128;

//Use the SmallVec implementation to save lots of memory 
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerNode{
    pub kmer: Kmer16,
    pub order: u32,
    pub order_val: u32,
    pub color: Color,
    pub child_nodes: SmallVec<[u32;1]>,
    pub child_edge_distance: SmallVec<[(u16,(Color,u8));1]>,
    pub parent_nodes: SmallVec<[u32;1]>,
    pub id: u32,
    pub canonical: bool,
    pub actual_ref_positions: SmallVec<[usize;0]>,
    pub repetitive: bool,
    pub primary_base: Option<u32>,
    pub closest_ref: u32,
    pub dist_to_closest_ref: u16,
}

pub struct Cmg<'a> {
    pub nodes: &'a mut Vec<KmerNode>,
    pub k: u8,
    pub properties: Vec<HashMap<String, u32>>,
}

impl<'a> Cmg<'a> {
    pub fn new(nodes: &mut Vec<KmerNode>, k: u8) -> Cmg {
        let mut new_cmg = Cmg {
            nodes,
            k,
            properties: Vec::new(),
        };
        new_cmg.properties.resize(new_cmg.nodes.len(), HashMap::new());
        new_cmg
    }

    pub fn out_degree(&self, node_id: u32) -> usize {
        self.nodes[node_id as usize].child_nodes.len()
    }

    pub fn in_degree(&self, node_id: u32) -> usize {
        self.nodes[node_id as usize].parent_nodes.len()
    }

    pub fn get_successors(&self, node_id: u32) -> Vec<u32> {
        self.nodes[node_id as usize].child_nodes.to_vec()
    }

    pub fn get_predecessors(&self, node_id: u32) -> Vec<u32> {
        self.nodes[node_id as usize].parent_nodes.to_vec()
    }

    pub fn get_property(&self, node_id: u32, property: String) -> Option<u32> {
        if let Some(value) = self.properties[node_id as usize].get(&property) {
            Some(*value)
        } else {
            None
        }
    }

    pub fn set_property(&mut self, node_id: u32, property: String, value: u32) {
        self.properties[node_id as usize].insert(property, value);
    }

    pub fn update_property<F>(&mut self, v: u32, property: String, func: F, values: Vec<u32>)
    where F: Fn(Vec<u32>) -> Option<u32> {
        if values.len() != 0 {
            if self.properties[v as usize].contains_key(&property) {
                let combined_vales = [self.properties[v as usize][&property]].iter().chain(values.iter()).cloned().collect::<Vec<u32>>();
                self.properties[v as usize].insert(property, func(combined_vales).unwrap());
            }
            else {
                if values.len() > 1 {
                    self.properties[v as usize].insert(property, func(values).unwrap());
                } else {
                    self.properties[v as usize].insert(property, values[0]);
                }
            }
        }
    }
    
}

pub struct ColoredCmg<'a> {
    pub graph: &'a mut Cmg<'a>,
}

impl ColoredCmg<'_> {
    pub fn new<'a>(graph: &'a mut Cmg<'a>) -> ColoredCmg<'a>
    where 'a : 'a {
        let new_colored_graph = ColoredCmg::<'a> {
            graph,
        };
        for node_id in 0..new_colored_graph.graph.nodes.len() {
            new_colored_graph.graph.set_property(node_id as u32, "color".to_string(), WHITE as u32);
        }
        new_colored_graph
    }

    pub fn set_color(&mut self, node_id: u32, color: u8) {
        // self.cmg.nodes[node_id as usize].properties.insert("color".to_string(), color as u32);
        self.graph.set_property(node_id, "color".to_string(), color as u32)
    }

    pub fn get_color(&self, node_id: u32) -> u8 {
        // self.cmg.nodes[node_id as usize].properties["color"] as u8
        self.graph.get_property(node_id, "color".to_string()).unwrap() as u8
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Bubble {
    pub id: u32,
    pub start: u32,
    pub end: u32,
    pub longest_path_length: Option<u32>,
    pub shortest_path_length: Option<u32>,
}

pub struct Graph {
    pub nodes: Vec<GraphNode>,
}

pub struct GraphNode {
    pub id: usize,
    pub children: Vec<usize>,
    pub parents: Vec<usize>,
    pub properties: HashMap<String, u32>,
}

impl Graph {
    pub fn out_degree(&self, node_id: usize) -> usize {
        self.nodes[node_id as usize].children.len()
    }

    pub fn in_degree(&self, node_id: usize) -> usize {
        self.nodes[node_id as usize].parents.len()
    }

    pub fn get_successors(&self, node_id: usize) -> &Vec<usize> {
        &self.nodes[node_id as usize].children
    }

    pub fn get_predecessors(&self, node_id: usize) -> &Vec<usize> {
        &self.nodes[node_id as usize].parents
    }

    pub fn get_property(&self, node_id: usize, property: String) -> Option<u32> {
        if let Some(value) = self.nodes[node_id as usize].properties.get(&property) {
            Some(*value)
        } else {
            None
        }
    }

    pub fn set_property(&mut self, node_id: usize, property: String, value: u32) {
        self.nodes[node_id as usize].properties.insert(property, value);
    }

    pub fn update_property<F>(&mut self, v: usize, property: String, func: F, values: Vec<u32>)
    where F: Fn(Vec<u32>) -> Option<u32> {
        let v_node = &mut self.nodes[v];
        if values.len() != 0 {
            if v_node.properties.contains_key(&property) {
                let combined_vales = [v_node.properties[&property]].iter().chain(values.iter()).cloned().collect::<Vec<u32>>();
                v_node.properties.insert(property, func(combined_vales).unwrap());
            }
            else {
                if values.len() > 1 {
                    v_node.properties.insert(property, func(values).unwrap());
                } else {
                    v_node.properties.insert(property, values[0]);
                }
            }
        }
    }
}

pub struct ColoredGraph<'a> {
    pub graph: &'a mut Graph
}

impl ColoredGraph<'_> {
    pub fn new(graph: &mut Graph) -> ColoredGraph {
        let mut new_colored_graph = ColoredGraph {
            graph
        };
        for node in new_colored_graph.graph.nodes.iter_mut() {
            node.properties.insert("color".to_string(), WHITE as u32);
        }
        new_colored_graph
    }

    pub fn set_color(&mut self, node_id: usize, color: u8) {
        self.graph.nodes[node_id as usize].properties.insert("color".to_string(), color as u32);
    }

    pub fn get_color(&self, node_id: usize) -> u8 {
        self.graph.nodes[node_id as usize].properties["color"] as u8
    }
}
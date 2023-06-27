use std::collections::HashMap;
use crate::constants::WHITE;

pub struct Graph {
    pub nodes: Vec<Node>,
}

pub struct Node {
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
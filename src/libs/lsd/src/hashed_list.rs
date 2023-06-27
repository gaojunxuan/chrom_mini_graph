use std::collections::HashMap;
use std::ops::{Index, Range, RangeInclusive};

/// A list implementation that allows for O(1) index_of lookup using a hash table
pub struct HashedList<T> {
    list: Vec<T>,
    dict: HashMap<T, usize>,
}

impl<T: Eq + std::hash::Hash + Clone> HashedList<T> {
    pub fn new(iterable: Option<impl IntoIterator<Item = T>>) -> Self {
        let mut hashed_list = HashedList {
            list: Vec::new(),
            dict: HashMap::new(),
        };
        
        if let Some(iter) = iterable {
            for v in iter {
                hashed_list.append(v);
            }
        }
        
        hashed_list
    }
    
    pub fn append(&mut self, v: T) {
        self.dict.insert(v.clone(), self.list.len());
        self.list.push(v);
    }
    
    pub fn index_of(&self, v: &T) -> isize {
        match self.dict.get(v) {
            Some(index) => *index as isize,
            None => -2,
        }
    }
    
    pub fn len(&self) -> usize {
        self.list.len()
    }
    
    pub fn get(&self, item: usize) -> Option<&T> {
        self.list.get(item)
    }
    
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.list.iter()
    }
    
    pub fn contains(&self, item: &T) -> bool {
        self.dict.contains_key(item)
    }
    
    pub fn add_pseudo_revers(&mut self, key: T, value: usize) {
        self.dict.insert(key, value);
    }
}

impl<T: std::fmt::Debug> std::fmt::Display for HashedList<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.list)
    }
}

impl<T: Eq + std::hash::Hash + Clone> Index<usize> for HashedList<T> {
    type Output = T;
    
    fn index(&self, index: usize) -> &Self::Output {
        &self.list[index]
    }
}

impl<T: Eq + std::hash::Hash + Clone> Index<Range<usize>> for HashedList<T> {
    type Output = [T];
    
    fn index(&self, index: Range<usize>) -> &Self::Output {
        &self.list[index]
    }
}

impl<T: Eq + std::hash::Hash + Clone> Index<RangeInclusive<usize>> for HashedList<T> {
    type Output = [T];
    
    fn index(&self, index: RangeInclusive<usize>) -> &Self::Output {
        &self.list[index]
    }
}

impl<T: Eq + std::hash::Hash + Clone> FromIterator<T> for HashedList<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut hashed_list = HashedList {
            list: Vec::new(),
            dict: HashMap::new(),
        };
        
        for v in iter {
            hashed_list.append(v);
        }
        
        hashed_list
    }
}
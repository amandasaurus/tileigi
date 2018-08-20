use std::rc::Rc;
use std::collections::HashSet;
use std::borrow::Borrow;
use std::ops::Deref;

pub struct StringStore(HashSet<Rc<String>>);

impl StringStore {
    pub fn new() -> Self {
        StringStore(HashSet::new())
    }

    pub fn get_string(&mut self, s: String) -> Rc<String> {
        let (result, should_add) = match self.0.get(&s) {
            None => {
                let rc: Rc<String> = Rc::new(s);
                (rc, true)
            },
            Some(value) => (value.clone(), false),
        };
        if should_add {
            self.0.insert(result.clone());
        }
        result
    }

    pub fn get_str(&mut self, s: &String) -> Rc<String> {
        let (result, should_add) = match self.0.get(s) {
            None => {
                let rc: Rc<String> = Rc::new(s.to_string());
                (rc, true)
            },
            Some(value) => (value.clone(), false),
        };
        if should_add {
            self.0.insert(result.clone());
        }
        result
    }
}

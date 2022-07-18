#[allow(non_snake_case)]
#[allow(unused_parens)]

mod nozzle;

#[no_mangle]
pub extern fn triple_input(input: i32) -> i32 {
    input * 3
}

// # References
// - https://www.engineeringtoolbox.com/methane-d_1420.html
// - https://en.wikipedia.org/wiki/File:Isentropic_Flow_Relations_Table.PNG
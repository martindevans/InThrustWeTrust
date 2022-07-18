/// Calculate the area ratio (area / throat_area)
///
/// # Arguments
/// - `mach`: mach number at exit
/// - `gamma`: ratio of specific heats of the gas
/// 
/// # Reference
/// - https://www.grc.nasa.gov/www/k-12/rocket/astar.html
/// - https://en.wikipedia.org/wiki/File:Isentropic_Flow_Relations_Table.PNG
#[no_mangle]
pub extern fn area_ratio_from_mach_number(mach: f64, gamma: f64) -> f64
{
    let gam = gamma;
    let M = mach;

    // A/A* = {[(gam+1)/2]^-[(gam+1)/(gam-1)/2]} / Me * [1 + Me^2 * (gam-1)/2]^[(gam+1)/(gam-1)/2] 

    // {[(gam+1)/2]^-[(gam+1)/(gam-1)/2]} / Me
    let left = f64::powf(
        (gam + 1.0) / 2.0,
        -((gam + 1.0) / (gam - 1.0) / 2.0)
    ) / M;

    // [1 + Me^2 * (gam-1)/2]^[(gam+1)/(gam-1)/2] 
    let right = f64::powf(
        (1.0 + f64::powf(M, 2.0) * (gam - 1.0) / 2.0),
        ((gam + 1.0)/(gam - 1.0) / 2.0)
    );

    return left * right;
}

/// Calculate the mass flow rate of the exhaust at the throat assuming choked flow
/// 
/// # Arguments
/// - `throat_area`: area of the nozzle at the throat
/// - `query_area`: area of the nozzle at the query location
/// - `gamma`: ratio of specific heats of the gas
/// - `r`: gas constant for gas
/// 
/// # Reference
/// - https://www.grc.nasa.gov/www/k-12/rocket/wcora.html
#[no_mangle]
pub extern fn nozzle_sonic_mass_flow_rate(throat_area: f64, pressure: f64, temperature: f64, gamma: f64, r: f64) -> f64
{
    let A = throat_area;
    let Pt = pressure;
    let Tt = temperature;
    let y = gamma;
    let R = r;

    let mdot = (Pt / f64::sqrt(R * Tt))
             * A
             * f64::sqrt(y)
             * f64::powf(
                (1.0 + (y - 1.0) / 2.0),
                (2.0 - 2.0 * y)
             );

    return mdot;
}

/// Calculate the total thrust for a given mass flow and exit velocity
#[no_mangle]
pub extern fn nozzle_thrust(mass_flow: f64, exit_velocity: f64) -> f64
{
    return mass_flow * exit_velocity;
}

#[no_mangle]
pub extern fn nozzle_exit_velocity(inlet_temperature: f64, inlet_pressure: f64, exit_pressure: f64, molecular_weight: f64, gamma: f64) -> f64
{
    // https://en.wikipedia.org/wiki/Rocket_engine_nozzle#de_Laval_nozzle_in_1_dimension

    let T = inlet_temperature;
    let R = 8314.5;
    let M = molecular_weight;
    let y = gamma;
    let pe = exit_pressure;
    let p = inlet_pressure;

    let pow_term = f64::powf(
        pe / p,
        (y - 1.0) / y
    );

    return f64::sqrt(
          ((T * R) / M)
        * ((2.0 * y) / (y - 1.0))
        * (1.0 - pow_term)
    );
}

#[cfg(test)]
mod tests
{
    #[test]
    fn nozzle_choked_flow_rate_raptor()
    {
        // Approximate numbers for the Raptor engine
        assert_eq!(
            f64::round(super::nozzle_sonic_mass_flow_rate(0.2216, 3000000.0, 3587.0, 1.4, 287.06)),
            670.0
        );
    }

    #[test]
    fn area_ratio_from_mach_number_solver()
    {
        // Numbers taken from solver on this page: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
        assert_eq!(
            f64::round(super::area_ratio_from_mach_number(2.940, 1.4)),
            4.0
        );
    }

    #[test]
    fn nozzle_exit_velocity_raptor()
    {
        // Numbers taken from example on this page: https://en.wikipedia.org/wiki/Rocket_engine_nozzle#de_Laval_nozzle_in_1_dimension
        assert_eq!(
            f64::round(super::nozzle_exit_velocity(
                3500.0,
                7000000.0,
                100000.0,
                22.0,
                1.22
            )),
            2802.0
        );
    }
}
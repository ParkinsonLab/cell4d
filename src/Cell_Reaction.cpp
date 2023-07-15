//=============================================================
// Cell_Reaction.cpp
// Donny Chan
#include "../inc/Cell_Reaction.h"
extern map<string, Species_Attributes *> species_details_map;
extern double timescale;
 

Cell_Reaction::Cell_Reaction() {
    
}

Cell_Reaction::~Cell_Reaction(){
    reactants_list.clear();
    products_list.clear();
    modifiers_list.clear();
}

// create a set of mappings that matches each reactant to the product. This assumes that no molecules are created or destroyed.
// the data type for both are int pairs that indicates index of the container, then index of the position within that container
void Cell_Reaction::set_reactant_product_pairing() {
    react_prod_map.clear();
    auto prod_list = products_list; // copy of product list to modify
    for(auto & react_container : reactants_list) {
        int reactant_index = 0;
        for(auto & reactant_mol : react_container.second) {
            bool react_mol_found = false;
            for(auto & prod_container : prod_list) {
                if(react_mol_found) break;
                int prod_index = 0;
                for(auto & prod_mol : prod_container.second) {
                    if(react_mol_found) break;
                    if(reactant_mol == prod_mol) {
                        react_prod_map[make_pair(react_container.first, reactant_index)] = make_pair(prod_container.first, prod_index);
                        prod_mol = "";
                        react_mol_found = true;
                    }
                    prod_index++;
                }
            }
            reactant_index++;
        }
    }
    cout << "matched reactants with products in reaction " << reaction_name << endl;
}

void Cell_Reaction::determine_interaction_type() { // by analyzing the reactant and product list of the reaction, determine how the reaction should proceed
    bool two_reactant_containers, two_product_containers, multi_product_containers;
    if(reactants_list.size() == 1) {
        two_reactant_containers = false;
    } else if(reactants_list.size() == 2) {
        two_reactant_containers = true;
    }

    // check product container numbers
    if(products_list.size() == 2) {
        two_product_containers = true, multi_product_containers = false;
    } else if(products_list.size() > 2) {
        two_product_containers = false, multi_product_containers = true;
    } else {
        two_product_containers = false, multi_product_containers = false;
    }


    if(two_reactant_containers) {

        if(two_product_containers) {
            bool no_species_change = false, matched_0 = false, matched_1 = false;
            if((reactants_list[0] == products_list[0] && reactants_list[1] == products_list[1]) || (reactants_list[1] == products_list[0] && reactants_list[0] == products_list[1])) {
                no_species_change = true, matched_0 = true, matched_1 = true;
            }
            if(!no_species_change) {
                if(reactants_list[0] == products_list[0] || reactants_list[0] == products_list[1]) matched_0 = true;
                if(reactants_list[1] == products_list[0] || reactants_list[1] == products_list[1]) matched_1 = true;
            }
            if(no_species_change) { // reactant and product species are the same
                interaction_type = "state_change"; // "cascade-style" reactions where only the states of the reactants are changed
            } else if(matched_0 != matched_1) { // one of the reacting containers don't match the reactants
                interaction_type = "modify"; // example A + B -> A + C
            } else { // both reactants are changed
                interaction_type = "new_products"; // example A + B -> C + D
            }
        } else if(multi_product_containers) { //
            interaction_type = "multi_product";
        } else if(!two_product_containers && !multi_product_containers) { // produces 1 container, could be attached to existing or a new one
            vector <string> total_reactants = reactants_list[0];
            total_reactants.insert(total_reactants.end(), reactants_list[1].begin(), reactants_list[1].end());
            if(is_permutation(total_reactants.begin(), total_reactants.end(), products_list[0].begin())) {
                interaction_type = "fusion";
            } else {
                interaction_type = "unknown";
            }
        }
    } else if(!two_reactant_containers) {
        if(two_product_containers) {
            interaction_type = "fission";
        } else if(multi_product_containers) {
            interaction_type = "split";
        } else {
            interaction_type = "state_change";
        }
    }

}

// determine interaction type for bulk_particle interactions
void Cell_Reaction::determine_bulk_interaction_type() {
    bool one_product = false;

    if(products_list.size() == 1) {
        one_product = true;
    }

    if(!one_product) {
        cerr << "bulk_particle interactions currently does not support reactions producing more than one product. exiting." << endl;
        exit(1);
    } else {
        interaction_type = "fusion";
    }

}

void Cell_Reaction::set_reaction_rate(bool AB_adjust) {
    if(forward_reaction_rate != 0) return;

    double difc_1 = 0, difc_2 = 0;
    for(auto & reactant : reactants_list[0]) {
        for(auto & species : species_details_map) {
            cout << species.first << endl;
        }
        cout << species_details_map[reactant]->diffusion_constant << endl;
        difc_1 += pow(species_details_map[reactant]->diffusion_constant, -3);
    }
    difc_1 = pow(difc_1, double(-1)/double(3));
    for(auto & reactant : reactants_list[1]) {
        difc_2 += pow(species_details_map[reactant]->diffusion_constant, -3);
    }
    difc_2 = pow(difc_2, double(-1)/double(3));

    if(AB_adjust) {
        // directly use fitted curves to get rate const from rstep (made of m_difc and binding radius)
        forward_reaction_rate = find_AB_rate_from_rstep((difc_1 + difc_2), binding_radius);
    } else {
        // apply standard smoluchowski's eq to get rate from binding radii and mutual difc
        forward_reaction_rate = AVO_CONST * 1000 * 4 * PI * (difc_1 + difc_2) * binding_radius;
    }

    cout << "Reaction " << reaction_name << " rate constant set to " << forward_reaction_rate << " from binding radius" << binding_radius << endl;
}


void Cell_Reaction::set_binding_radius(bool AB_adjust) {
    if(binding_radius != 0) return;

    double difc_1 = 0, difc_2 = 0;
    for(auto & reactant : reactants_list[0]) {
        difc_1 += pow(species_details_map[reactant]->diffusion_constant, -3);
    }
    difc_1 = pow(difc_1, double(-1)/double(3));
    for(auto & reactant : reactants_list[1]) {
        difc_2 += pow(species_details_map[reactant]->diffusion_constant, -3);
    }
    difc_2 = pow(difc_2, double(-1)/double(3));

    double mutual_difc = difc_1 + difc_2;

    if(AB_adjust) {
        // AB-adjusted Smoluchowski equation using two fitted curves
        binding_radius = calculate_AB_radius(forward_reaction_rate, mutual_difc);
    } else {
        // Smoluchowski's equation here
        // forward_reaction_rate = 4 * PI * (difc_1 + difc_2) * binding_radius;
        binding_radius = forward_reaction_rate / (AVO_CONST * 1000 * 4 * PI * mutual_difc);
    }
    
    cout << "Reaction " << reaction_name << " binding radius set to " << binding_radius << " from rate constant" << forward_reaction_rate << endl;

}

// for reactions with more than one product
void Cell_Reaction::set_unbinding_radius(double unbinding_radius) {
    unbind_rad_min = unbinding_radius;
    unbind_rad_max = unbinding_radius * 1.2;
    cout << "Reaction " << reaction_name << " unbinding radius set to " << unbind_rad_min << endl;
}

// using a fitted curve based on AB-adjusted Smoluchowski radius, calculate the correct binding radius for
// reactants given their mutual diffusion constant and the target reaction rate constant
double Cell_Reaction::calculate_AB_radius(double target_r_rate, double m_difc) {
    double rate_diff = 0;
    double radius_guess = target_r_rate / (AVO_CONST * 1000 * 4 * PI * m_difc); // Inverted Smoluchowski's eq for radius
    double guess_rate = 0;
    bool found_rate = false;
    int counter = 0;
    double adjustment = 0;

    // iteratively try to find the correct radius for desired reaction rate constant
    while(!found_rate) {
        guess_rate = find_AB_rate_from_rstep(m_difc, radius_guess);
        rate_diff = (guess_rate - target_r_rate) / target_r_rate;
        adjustment = rate_diff / 2;

        // check if guess is close enough yet
        if(abs(rate_diff) > 0.000005) {
            // if guess is too small, increase by half the error between current and target
            // if guess is too big, divide current guess by error
			if(rate_diff < 1) {
				radius_guess = radius_guess + (radius_guess * -adjustment);
			} else {
				radius_guess = radius_guess / rate_diff;
			}
        } else {
            // if guessed radius gives a rate that's within threshold of 0.000005, stop searching
            found_rate = true;
        }
        counter += 1;
        if(counter > 100) {
            cerr << "Warning: radius not found after 100 tries, will proceed with current value" << endl;
            cerr << "Guess rate: " << guess_rate << " intended: " << target_r_rate << " diff:" << rate_diff << endl;
            found_rate = true;
        }
    }

    // final check on the reduced step length. If it is less than the smallest value of lookup table,
    // adjustment is not necessary. reset to smoluchowski radius.
    double final_reduced_step = sqrt(2 * m_difc * timescale) / radius_guess;
    if(final_reduced_step < 0.0498) {
        radius_guess = target_r_rate / (AVO_CONST * 1000 * 4 * PI * m_difc);
    }

    return radius_guess;
}

// using fitted 5-parameter logistic function to AB-adjusted smoluchowski table, find rate
double Cell_Reaction::find_AB_rate_from_rstep(double m_difc, double bind_radius) {
    double mutual_reduced_step = sqrt(2 * m_difc * timescale) / bind_radius;

    // for steps smaller than the minimum of the lookup table, the curve is not valid.
    // since small reduced steps converge to smoluchowski's rate anyways, just return that.
    if(mutual_reduced_step < 0.0498) {
        double rate_const = 4 * PI * m_difc * bind_radius * AVO_CONST * 1000;
        return(rate_const);
    }

    double b, c, d, e, f;
    // 2 curves fitted to the lookup table. Use one set for r_step > 1, other one for r_step < 1
    if(mutual_reduced_step > 1) {
		b = -2.89171179;
		c = -0.00648886;
		d = 4.19026581;
		e = 0.00279511;
		f = 0.62295397;
    } else {
		b = -2.35948390;
		c = -0.00125267;
		d = 4.78882764;
		e = 0.00154475;
		f = 0.80220128;
    }
    double exponent_base_term = 1 + exp(b * (log(mutual_reduced_step) - e));
    double reduced_rate = c + (d - c) / (pow(exponent_base_term, f));

    double rate_const = reduced_rate * pow(bind_radius, 3) * AVO_CONST * 1000 / timescale;
    return(rate_const);
}


string Cell_Reaction::get_interaction_type() {
    return interaction_type;
}


bool Cell_Reaction::is_two_reactants() {
    if(reactants_list.size() == 2) return true;
    return false;
}

bool Cell_Reaction::is_fission() {
    if(products_list.size() > 1) return true;
    return false;
}

function X_final=equivalance_ratio(global_stoic,fuel_matrix,oxidizer_matrix,phi,oxidizer_divide,fuel_divide)
total_fuel=phi.*fuel_matrix./str2double(fuel_divide{2});
total_oxidizer=oxidizer_matrix./str2double(oxidizer_divide{2});
X_fuel=sum(total_fuel.*global_stoic)/(sum(total_fuel.*global_stoic)...
    +sum(total_oxidizer.*global_stoic));
X_air=1-X_fuel;
final_X_fuel=X_fuel.*fuel_matrix;
final_X_air=X_air.*oxidizer_matrix;
X_final=final_X_fuel+final_X_air;
end
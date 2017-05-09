extern void initialize();
extern void initialize_files();
extern void initialize_command_line_variables();

extern void symmetric_initial_condition1();
extern void initial_condition1();
extern void initial_condition2();
extern void initial_condition3();
extern void initial_condition4();
extern void standard_kt_initial_condition();
extern void symmetric_stretched_kt_initial_condition();
extern void stretched_kt_initial_condition();

extern void evolve_mts();
extern void calc_rates();

extern void evolve_kts(unsigned long long);
extern void mt_kt_int();
extern void dash_pot(double, double, double);

extern void update_system();
extern void set_prev_mono_conc(double);
extern void check_kinetochore_position();
extern void shift_simulation(bool);


extern void print_mt_data();
extern void print_mean_mt_length();
extern void print_tension();
extern void print_distrib_stats();

extern void print_kt_data();
extern void print_kt_pos();

extern void write_restart();
extern void read_restart();
extern void perturb();


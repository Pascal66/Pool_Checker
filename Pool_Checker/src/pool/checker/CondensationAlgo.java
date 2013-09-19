package pool.checker;

import android.util.Log;

public class CondensationAlgo {
    // http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/ISARD1/condensation.html
    // "model_parameters.h"
    /** The following are all the constants controlling the behaviour of the system and output. */
    /** How many samples in the distribution? */
    int NSamples = 1000;

    /** How many iterations to run the filter? */
    int NIterations = 100;

    /** The simulated object follows a model of the same form as the process */
    double SimulatedMean = -0.1;
    double SimulatedScaling = 0.4;
    double SimulatedSigma = 0.075;
    double SimulatedMeasSigma = 0.03;

    /** The prior distribution over samples at the first timestep is a Gaussian with the following parameters. */
    double PriorMean = 0.0;
    double PriorSigma = 0.2;

    /** The process model is a first-order Auto-Regressive Process of the following form:
     * x_{t+1} - ProcessMean =
     * (x_t - ProcessMean) * ProcessScaling + ProcessSigma * w_t
     * where w_t is zero-mean unit iid Gaussian noise. */
    double ProcessMean = -0.1;
    double ProcessScaling = 0.4;
    double ProcessSigma = 0.075;

    /** The observation density is a mixture of Gaussians, where each observed object has a different sigma as follows. */
    double ObsSigma = 0.03;

    /** How many columns wide is the ASCII histogram? */
    int HistBins = 79;
    /** How many lines of text does the ASCII histogram take? */
    int HistLines = 25;
    /**  What is the highest value the density histogram can represent before saturating above HistLines?  */
    double MaxHistHeight = 0.2;

    /** What is the distance from the origin to the edge of the histogram? */
    double DisplayWidth = 0.35;
    
    // "data_types.h"
    /* The following are model-specific definitions */

    /** The state vector for this simple model is one-dimensional. In multiple dimensions, for example, Sample would be an array.
     */
    double[] Sample;

    /** The process model used in this simple implementation is a one-dimensional first-order auto-regressive process, where
     * 
     * (x_t - mean) = (x_{t-1} - mean) * scaling + sigma w_t
     * 
     * and the w_t are iid unit zero-mean Gaussian noise samples.
     */
    class ProcessModel {
	double mean;
	double scaling;
	double sigma;
    }

    /**
     * This structure contains the parameters of the simulation model. The
     * simulated data is produced from a model of a single particle
     * following a one-dimensional 1st-order ARP whose parameters are
     * stored in the structure process. Measurements are simulated by
     * adding Gaussian noise with std. deviation sigma to the true
     * simulated position of the particle.
     */
    class SceneModel {
	ProcessModel process = new ProcessModel();
	double sigma;
    }

    /**
     * This structure contains the data for the measurements made at each
     * timestep. The simulation consists of a single point-measurement.
     * This structure stores the true position of the simulated particle,
     * and its measured position, which is corrupted by noise.
     */
    class MeasData {
	double truePosition, observed;
    }

    /**
     * The prior distribution of the state is taken to be Gaussian with
     * the parameters stored in this structure.
     */
    class PriorModel {
	double mean, sigma;
    }

    /**
     * The observation model is of Gaussian noise with std. deviation
     * sigma, so
     * 
     * z_t = x_t + sigma v_t
     * 
     * where the v_t are assumed iid unit zero-mean Gaussian noise
     * samples. In principle for this simulation, the modelled observation
     * noise can be different to the simulated observation noise, hence
     * the presence of a separate std. deviation parameter in the
     * structure SceneModel.
     */
    class ObservationModel {
	double sigma;
    }

    /**
     * This structure contains information about how the state
     * distribution should be displayed at each timestep. The display in
     * this implementation is a one-dimensional histogram, and this
     * parameter sets the width of the histogram (so that the interval
     * [-histogram_width, histogram_width] is displayed).
     */
    class DisplayModel {
	double histogram_width;
    }

    /* End of model-specific structures */

    /* The following should work with any model-specific definitions */

    /**
     * This structure contains all the parameter settings for the models
     * which remain constant throughout a run of the algorithm.
     */
    class GlobalData {
	/** The parameters specifying the model of the prior distribution for the first timestep. */
	PriorModel prior = new PriorModel();
	/** The parameters specifying the process model. */
	ProcessModel process = new ProcessModel();
	/** The parameters specifying the simulation model of the scene. This is only used in the case of simulated data. */
	SceneModel scene = new SceneModel();
	/** The parameters specifying the observation model. */
	ObservationModel obs = new ObservationModel();
	/** The parameters specifying how to display the state estimate at each timestep. */
	DisplayModel disp = new DisplayModel();
	/** The number of samples N. */
	int nsamples;
	/** The number of iterations to run the filter. */
	int niterations;
    }

    /**
     * This structure contains all of the information which is specific to
     * a given iteration of the algorithm.
     */
    class IterationData {
	double[] new_positions;

	/**
	 * The following arrays contain the sample positions for the current
	 * and previous timesteps respectively. At the end of each
	 * iteration, these pointers are swapped over to avoid copying data
	 * structures, so their addresses should not be relied on.
	 */
	double[] old_positions;

	/**
	 * The following arrays give the sample weights and cumulative
	 * probabilities, as well as the largest cumulative
	 * probability. There is no stage in the algorithm when the weights
	 * from both the previous and current timesteps are needed. At the
	 * beginning of an iteration, sample_weights contains the weights
	 * from the previous iteration, and by the end it contains the
	 * weights of the current iteration. The cumulative probabilities
	 * are not normalised, so largest_cumulative_prob is needed to store
	 * the largest cumulative probability (for simplicity of the binary
	 * search algorithm, cumul_prob_array[0] = 0)
	 */
	double[] sample_weights, cumul_prob_array;
	double largest_cumulative_prob;

	/**
	 * The measurements made in a given iteration are stored here. For
	 * some applications a discrete set of measurements is not
	 * appropriate, and this could contain, e.g. a pointer to an image
	 * structure.
	 */
	MeasData meas = new MeasData();
    }

    /* End of generic structures */
    /* End of global variables */

    // #include "condensation.h"

    /* All of the global information is packaged into the following two
     * structures. `global' contains information which is constant over a
     * run of the algorithm and `data' contains all of the current state
     * at a given iteration. */

    GlobalData global = new GlobalData();
    IterationData data = new IterationData();

    /* End of global variables */

    /* From here on is generic Condensation regardless of the form of
     * model or observation used. All of the model-specific routines are
     * found in model_specific.c */

    /**
     * This is binary search using cumulative probabilities to pick a base
     * sample. The use of this routine makes Condensation O(NlogN) where N
     * is the number of samples. It is probably better to pick base
     * samples deterministically, since then the algorithm is O(N) and
     * probably marginally more efficient, but this routine is kept here
     * for conceptual simplicity and because it maps better to the
     * published literature.
     */
    int pick_base_sample() {
	double choice = uniform_random() * data.largest_cumulative_prob;
	int low, middle, high;

	low = 0;
	high = global.nsamples;

	while (high > (low + 1)) {
	    middle = (high + low) / 2;
	    if (choice > data.cumul_prob_array[middle])
		low = middle;
	    else
		high = middle;
	}

	return low;
    }

    /**
     * This routine computes all of the new (unweighted) sample
     * positions. For each sample, first a base is chosen, then the new
     * sample position is computed by sampling from the prediction density
     * p(x_t|x_t-1 = base). predict_sample_position is obviously
     * model-dependent and is found in model_specific.c, but it can be
     * replaced by any process model required.
     */
    void predict_new_bases() {
	int n, base;

	for (n = 0; n < global.nsamples; ++n) {
	    base = pick_base_sample();
	    predict_sample_position(n, base);
	}
    }

    /**
     * Once all the unweighted sample positions have been computed using
     * predict_new_bases, this routine computes the weights by evaluating
     * the observation density at each of the positions. Cumulative
     * probabilities are also computed at the same time, to permit an
     * efficient implementation of pick_base_sample using binary
     * search. evaluate_observation_density is obviously model-dependent
     * and is found in model_specific.c, but it can be replaced by any
     * observation model required.
     */
    void calculate_base_weights() {
	int n;
	double cumul_total;

	cumul_total = 0.0;
	for (n = 0; n < global.nsamples; ++n) {
	    data.sample_weights[n] = evaluate_observation_density(n);
	    data.cumul_prob_array[n] = cumul_total;
	    cumul_total += data.sample_weights[n];
	}
	data.largest_cumulative_prob = cumul_total;
    }

    /**
     * Go and output the estimate for this iteration (which is a
     * model-dependent routine found in model_specific.c) and then swap
     * over the arrays ready for the next iteration.
     */
    void update_after_iterating(int iteration) {
	double[] /* Sample */temp;

	display_data(iteration);

	temp = data.new_positions;
	data.new_positions = data.old_positions;
	data.old_positions = temp;
    }

    /**
     * / obtain_observations is model-dependent and can be found in model_specific.c
     */
    void run_filter() {
	int i;

	for (i = 0; i < global.niterations; ++i) {
	    obtain_observations(); /* Go make necessary measurements */
	    predict_new_bases(); /* Push previous state through process model */
	    calculate_base_weights(); /* Apply Bayesian measurement weighting */
	    update_after_iterating(i); /* Tidy up, display output, etc. */
	}
    }

    /**
     * This routine fills in the data structures with default constant
     * values. It could be enhanced by reading informatino from the
     * command line to allow e.g. N to be altered without recompiling.
     */
    void initialise_defaults() {
	global.nsamples = NSamples;
	global.niterations = NIterations;

	initialise_model_specific_defaults();
    }

    /**
     * Create all the arrays, then fill in the prior distribution for the
     * first iteration. The prior is model-dependent, so
     * set_up_prior_conditions can be found in model_specific.c
     */
    boolean initialise() {
	initialise_defaults();

	data.new_positions = new double[(/* sizeof(Sample) * */global.nsamples)];
	data.old_positions = new double[(/* sizeof(Sample) * */global.nsamples)];

	data.sample_weights = new double[(/* sizeof(double) * */global.nsamples)];
	data.cumul_prob_array = new double[(/* sizeof(double) * */global.nsamples)];

	// if (!data.new_positions || !data.old_positions || !data.sample_weights || !data.cumul_prob_array) {
	// fprintf(stderr, "Failed to allocate memory for sample arrays\n");
	// return false;
	// }

	set_up_prior_conditions();

	return true;
    }

    /* Tidy up */
    void shut_down() {
	// free(data.new_positions);
	// free(data.old_positions);
	// free(data.sample_weights);
	// free(data.cumul_prob_array);
    }

    boolean main() {
	if (!initialise())
	    return true;

	run_filter();

	shut_down();

	return false;
    }

    // #include "data_types.h"
    // #include "model_parameters.h"
    // #include "condensation.h"

    /* The following routines are model-specific and should be replaced to
     * implement an arbitrary process and observation model. */

    void initialise_model_specific_defaults() {
	/* Set up the parameters of the simulation model, process and
	 * observation. */
	global.scene.process.mean = SimulatedMean;
	global.scene.process.scaling = SimulatedScaling;
	global.scene.process.sigma = SimulatedSigma;
	global.scene.sigma = SimulatedMeasSigma;

	/* Set up the parameters of the prior distribution */
	global.prior.mean = PriorMean;
	global.prior.sigma = PriorSigma;

	/* Set up the parameters of the process model */
	global.process.mean = ProcessMean;
	global.process.scaling = ProcessScaling;
	global.process.sigma = ProcessSigma;

	/* Set up the parameters of the observation model */
	global.obs.sigma = ObsSigma;

	/* Set up the parameters of the display model */
	global.disp.histogram_width = DisplayWidth;
    }

    /**
     * Set up the initial sample-set according to the prior model. The
     * prior is a simple Gaussian, so each of the positions is filled in
     * by sampling that Gaussian and the weights are initialised to be
     * uniform. gaussian_random can be found in utility.c
     */
    void set_up_prior_conditions() {
	int n;

	for (n = 0; n < global.nsamples; ++n) {
	    data.old_positions[n] = global.prior.mean + global.prior.sigma
		    * gaussian_random();

	    /* The probabilities are not normalised. */
	    data.cumul_prob_array[n] = (double) n;
	    data.sample_weights[n] = 1.0;
	}

	/**
	 * The probabilities are not normalised, so store the largest value
	 * here (for simplicity of the binary search algorithm,
	 * cumul_prob_array[0] = 0). This can then be used as a
	 * multiplicative normalisation constant for the sample_weights
	 * array, as well.
	 */
	data.largest_cumulative_prob = (double) n;

	/* This is the initial positions the simulated object. */
	data.meas.truePosition = 0.0;
    }

    /**
     * The process model for a first-order auto-regressive process is:
     * 
     * x_{t+1} - mean = (x_t - mean)*scaling + sigma*w_t
     * 
     * where w_t is unit iid Gaussian noise. gaussian_random can be found
     * in utility.c
     */
    double iterate_first_order_arp(double truePosition, ProcessModel process) {
	return process.mean + ((truePosition - process.mean) * process.scaling)
		+ process.sigma * gaussian_random();
    }

    /**
     * In a real implementation, this routine would go and actually make
     * measurements and store them in the data.meas structure. This
     * simulation consists of an `object' moving around obeying a
     * first-order auto-regressive process, and being observed with its
     * true positions corrupted by Gaussian measurement noise.
     * Accordingly, this routine calculates the new simulated true and
     * measured position of the object.
     */
    void obtain_observations() {
	data.meas.truePosition = iterate_first_order_arp(
		data.meas.truePosition, global.scene.process);

	data.meas.observed = data.meas.truePosition + global.scene.sigma
		* gaussian_random();
    }

    /**
     * This routine samples from the distribution
     * 
     * p(x_t | x_{t-1} = old_positions[old_sample])
     * 
     * and stores the result in new_positions[new_sample]. This is
     * straightforward for the simple first-order auto-regressive process
     * model used here, but any model could be substituted.
     */
    void predict_sample_position(int new_sample, int old_sample) {
	data.new_positions[new_sample] = iterate_first_order_arp(
		data.old_positions[old_sample], global.process);
    }

    /**
     * This routine evaluates the observation density
     * 
     * p(z_t|x_t = new_positions[new_sample])
     * 
     * The observation model in this implementation is a simple mixture of
     * Gaussians, where each simulated object is observed as a 1d position
     * and measurement noise is represented as Gaussian. For a
     * visual-tracking application, this routine would go and evaluate the
     * likelihood that the object is present in the image at the position
     * encoded by new_positions[new_sample]. evaluate_gaussian can be
     * found in utility.c
     */
    double evaluate_observation_density(int new_sample) {
	return evaluate_gaussian(data.new_positions[new_sample]
		- data.meas.observed, global.obs.sigma);
    }

    // /**
    // * The following two routines provide rudimentary graphical output in
    // * the form of an ASCII histogram of the 1d state density.
    // */
    // int eval_bin(float where) {
    // return (HistBins - 1)
    // / 2
    // + (int) (where * (HistBins - 1) / 2 / global.disp.histogram_width);
    // }

    // void display_histogram(double estimated_position)
    // {
    // double bins[HistBins];
    //
    // int b, n, line, which_bin, meas_bin, est_bin, true_bin;
    // double lineheight;
    // char outc;
    //
    // for (b=0; b<HistBins; ++b)
    // bins[b] = 0.0;
    //
    // for (n=0; n<global.nsamples; ++n) {
    // which_bin = eval_bin(data.new_positions[n]);
    //
    // if (which_bin >=0 && which_bin < HistBins)
    // bins[which_bin] +=
    // data.sample_weights[n] / data.largest_cumulative_prob;
    // }
    //
    // for (b=0; b<HistBins; ++b)
    // bins[b] = (bins[b] * (double) HistLines) / MaxHistHeight;
    //
    // for (line=0; line<HistLines; ++line) {
    // lineheight = (double) (HistLines-1-line);
    // for (b=0; b<HistBins; ++b) {
    // if (line==0 && bins[b] >= lineheight+1.0)
    // outc = '*';
    // else if (bins[b] >= lineheight+0.5 && bins[b] < lineheight+1.0)
    // outc = '-';
    // else if (bins[b] >= lineheight && bins[b] < lineheight+0.5)
    // outc = '_';
    // else
    // outc = ' ';
    // printf("%c", outc);
    // }
    // printf("\n");
    // }
    //
    // true_bin = eval_bin(data.meas.true);
    // meas_bin = eval_bin(data.meas.observed);
    // est_bin = eval_bin(estimated_position);
    //
    // for (b=0; b<HistBins; ++b) {
    // if ((b == meas_bin && b == est_bin) ||
    // (b == meas_bin && b == true_bin) ||
    // (b == true_bin && b == est_bin))
    // outc = '*';
    // else if (b == true_bin)
    // outc = '.';
    // else if (b == meas_bin)
    // outc = '+';
    // else if (b == est_bin)
    // outc = 'x';
    // else
    // outc = ' ';
    // printf("%c", outc);
    // }
    // printf("\n");
    // }

    /**
     * This routine computes the estimated position of the object by
     * estimating the mean of the state-distribution as a weighted mean of
     * the sample positions, then displays a histogram of the state
     * distribution along with a character denoting the estimated position
     * of the object.
     */
    void display_data(int iteration) {
	int n;
	double aggregate;

	aggregate = 0.0;

	/* Compute the unnormalised weighted mean of the sample
	 * positions. */
	for (n = 0; n < global.nsamples; ++n)
	    aggregate += data.new_positions[n] * data.sample_weights[n];

	aggregate /= data.largest_cumulative_prob;

	// display_histogram(aggregate);

	Log.i("Condensation", iteration + ": Measured pos. "
		+ data.meas.observed + " TruePos " + data.meas.truePosition
		+ " Est. position " + aggregate);

	// #ifdef ANSI_TERM_SEQUENCES
	// /* If possible, run the cursor back up the screen so the histogram
	// stays in the same place instead of scrolling down the display. */
	// printf("\033[%dA", HistLines+2);
	// sleep(1);
	// #endif
    }

    /* End of model-specific routines */

    /* The following are some utility routines. */

    /**
     * This is the worst random number routine in the known universe, but
     * I use it for portability. Feel free to replace it.
     */
    double uniform_random() {
	return Math.random()/* / RAND_MAX */;
    }

    /**
     * This Gaussian routine is stolen from Numerical Recipes and is their
     * copyright.
     */

    double gaussian_random() {
	int next_gaussian = 0;
	double saved_gaussian_value = 0;

	double fac, rsq, v1, v2;

	if (next_gaussian == 0) {
	    do {
		v1 = 2.0 * uniform_random() - 1.0;
		v2 = 2.0 * uniform_random() - 1.0;
		rsq = v1 * v1 + v2 * v2;
	    } while (rsq >= 1.0 || rsq == 0.0);
	    fac = Math.sqrt(-2.0 * Math.log(rsq) / rsq);
	    saved_gaussian_value = v1 * fac;
	    next_gaussian = 1;
	    return v2 * fac;
	} else {
	    next_gaussian = 0;
	    return saved_gaussian_value;
	}
    }

    double evaluate_gaussian(double val, double sigma) {
	/** This private definition is used for portability */
	double Condense_PI = 3.14159265358979323846;

	return 1.0 / (Math.sqrt(2.0 * Condense_PI) * sigma)
		* Math.exp(-0.5 * (val * val / (sigma * sigma)));
    }

    /* End of utility routines */

}

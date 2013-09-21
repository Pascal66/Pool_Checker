package net.priveyes.ia.controller;

import net.priveyes.fuzzy.engine.FuzzyBlockOfRules;
import net.priveyes.fuzzy.engine.FuzzyEngine;
import net.priveyes.fuzzy.engine.LinguisticVariable;
import net.priveyes.fuzzy.engine.NoRulesFiredException;

public class FuzzyController {
	private Double result;

	private FuzzyEngine fuzzyEngine = null;
	private FuzzyBlockOfRules fuzzyRules = null;

	LinguisticVariable angle = null;
	LinguisticVariable turn = null;

	// Tactic Variables
	private String[] rules;

	public Double fuzzyControl(Double steering) {
//		fuzzyEngine.reset();
		// 1. Create Lingustic variables and define membership functions
		angle = new LinguisticVariable("angle");
		angle.add("negative", -314, -3.14, -1, 0);
		angle.add("positive", 0, 1, 3.14, 314);
		turn = new LinguisticVariable("turn");
		turn.add("left", -2, -1, -1, 0);
		turn.add("right", 0, 1, 1, 2);
		// 2. Create a fuzzy engine
		fuzzyEngine = new FuzzyEngine();
		// 3. Register all LVs
		fuzzyEngine.register(angle);
		fuzzyEngine.register(turn);

		// 4. Create a block of rules
		rules = new String[] { "if angle is negative then turn is left",
				"if angle is positive then turn is right" };
		fuzzyRules = new FuzzyBlockOfRules(rules);

		// 5. Register the block
		fuzzyEngine.register(fuzzyRules);

		// Feeding input values to the Linguistic Variables
		angle.setInputValue(steering);

		// 6. Parse the rules
		try {
			fuzzyRules.parseBlock();
		} catch (Exception e) {
			System.out.print(e.getMessage() + "\n");
		}
		// 7. Perform the evaluation
		try {
			fuzzyRules.evaluateBlock();
			// - slower execution, returns a String with evaluation results for
			// every fuzzy expression
			// Log.w(TAG, fuzzyBlockOfRules.evaluateBlockText());
		} catch (Exception e) {
			System.out.print(e.getMessage() + "\n");
		}
		// 8. Obtain the result(s)
		try {
			result = turn.defuzzify();
		} catch (NoRulesFiredException nrfe) {
			return 0.0;
		}
		return result;

	}
}
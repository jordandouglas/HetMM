package hetmm;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;


@Description("A beta distribution upper and lower bound between the product of two variances")
public class CovarianceBetaPrior extends Distribution {
	
	
	final public Input<RealParameter> alphaInput = new Input<>("alpha", "beta distribution alpha parameter", Input.Validate.REQUIRED);
	final public Input<RealParameter> betaInput = new Input<>("beta", "beta distribution beta parameter", Input.Validate.REQUIRED);
	
	final public Input<RealParameter> mean1Input = new Input<>("mean1", "mean of parameter 1", Input.Validate.REQUIRED);
	final public Input<RealParameter> mean2Input = new Input<>("mean2", "mean of parameter 2", Input.Validate.REQUIRED);
	
	final public Input<RealParameter> rhoInput = new Input<>("rho", "correlation between two parameters", Input.Validate.REQUIRED);
	final public Input<RealParameter> var1Input = new Input<>("var1", "variance of paramreter 1", Input.Validate.REQUIRED);
	final public Input<RealParameter> var2Input = new Input<>("var2", "variance of parameter 2", Input.Validate.REQUIRED);

	
	public double calculateLogP() {
        logP = 0;
        
        double alpha = alphaInput.get().getValue();
        double beta = betaInput.get().getValue();
        double mean1 = mean1Input.get().getValue();
        double mean2 = mean2Input.get().getValue();
        double var1 = var1Input.get().getValue();
        double var2 = var2Input.get().getValue();
        double rho = rhoInput.get().getValue();
        
        double covar = 	rho*Math.sqrt(var1)*Math.sqrt(var2);
        
        // What is min and max of covar?
        double limit = var1*var2;
        
        if (Math.abs(covar) > limit) {
        	logP = Double.NEGATIVE_INFINITY;
        	return logP;
        }
        
        
        // Normalise covar into this range
        double x = (covar + limit) / (2*limit);
        
        
        // Confirm that the MVN is valid
 		double[] means = new double[] { Math.log(mean1), Math.log(mean2)};
 		double[][] variances = new double[2][2];
 		variances[0][0] = var1;
 		variances[0][1] = covar;
 		variances[1][0] = covar;
 		variances[1][1] = var2;
 		try {
 			MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(means, variances);
 		}catch (Exception e) {
 			logP = Double.NEGATIVE_INFINITY;
        	return logP;
 		}
        
        BetaDistribution dist = new BetaDistribution(alpha, beta);
        logP = dist.logDensity(x);
        return logP;
    }
	
	
	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(rhoInput.get().getID());
		return args;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(var1Input.get().getID());
		conds.add(var2Input.get().getID());
		conds.add(alphaInput.get().getID());
		conds.add(betaInput.get().getID());
		return conds;
	}

	@Override
	public void sample(State state, Random random) {
		
		if (this.sampledFlag) return;
		this.sampleConditions(state, random);
		
		double alpha = alphaInput.get().getValue();
        double beta = betaInput.get().getValue();
        double var1 = var1Input.get().getValue();
        double var2 = var2Input.get().getValue();
        double mean1 = mean1Input.get().getValue();
        double mean2 = mean2Input.get().getValue();
	        
        BetaDistribution dist = new BetaDistribution(alpha, beta);
		
        
        // Try again if the MVN is invalid for whatever reason
        int nattempts = 0;
        while (true) {
        
			// Sample from beta
			double x = dist.sample(); // TODO use random 
			
			// Normalise
			double limit = Math.sqrt(var1)* Math.sqrt(var2);
			double covar = x*2*limit - limit;
			
			double rho = covar/limit;
			
			
			// Confirm that the MVN is valid
			double[] means = new double[] { Math.log(mean1), Math.log(mean2)};
	 		double[][] variances = new double[2][2];
	 		variances[0][0] = var1;
	 		variances[0][1] = covar;
	 		variances[1][0] = covar;
	 		variances[1][1] = var2;
	 		boolean valid = true;
	 		try {
	 			MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(means, variances);
	 		}catch (Exception e) {
	 			valid = false;
	 			
	 			if (nattempts > 10000) {
	 				e.printStackTrace();
		 			throw new IllegalArgumentException("Unable to sample covariance because the MVN is still throwing exceptions after 10,000 attemps");
		 		}
	 		}
	 		
	 		nattempts ++;
	 		
	 		if (valid) {
	 			
	 			// Set value
	 			rhoInput.get().setValue(rho);
	 			break;
	 		}
	 		
		
        }
		
		
		
		this.sampledFlag = true;
		
	}

}

package hetmm;

import java.io.PrintStream;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Likelihood for fitting a single enzyme to the Michaelis-Menten model")
public class HomogeneousLikelihood extends Distribution {
	
	final public Input<KineticData> dataInput = new Input<>("data", "the data", Input.Validate.OPTIONAL);
	
	final public Input<RealParameter> kmInput = new Input<>("km", "michaelis constant Km", Input.Validate.REQUIRED);
	final public Input<RealParameter> vmaxInput = new Input<>("vmax", "limiting rate Vmax", Input.Validate.REQUIRED);
	final public Input<RealParameter> sdInput = new Input<>("sd", "standard deviation of error", Input.Validate.REQUIRED);
	
	
	// Hill coefficient?
	final public Input<RealParameter> hillInput = new Input<>("hill", "hill lower coefficient", Input.Validate.OPTIONAL);
	final public Input<IntegerParameter> hillIndicatorInput = new Input<>("hillIndicator", "hill coefficient indicator. 0 is h=1, -1 is h<1, +1 is h>1", Input.Validate.OPTIONAL);
	
	final public Input<Boolean> logInput = new Input<>("log", "perform regression in log space?", true);
	
	
	
	@Override
	public void initAndValidate() {
		
		if (hillIndicatorInput.get() != null) {
			
			if (hillInput.get() == null) {
				throw new IllegalArgumentException("Please provide both hill inputs, or none at all. hillL and hillIndicator");
			}
			
			
			hillIndicatorInput.get().setLower(-1);
			hillIndicatorInput.get().setUpper(+1);
			
		}
		
		

	}
	
	
	@Override
	public double calculateLogP() {
		
		logP = 0;
		
		KineticData data = dataInput.get();
		if (data == null) {
			return logP;
		}
		
		boolean valid = this.prepare();
		if (!valid) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}
		double error;
		double sd = sdInput.get().getValue();
		NormalDistribution dist = new NormalDistributionImpl(0, sd);
		
		// Iterate through observations
		for (int i = 0; i < data.getNumberOfObservations(); i ++) {
			
			
			// Is there data?
			if (!data.hasV(i)) {
				continue;
			}
			
			// Observations
			double vObs = data.getV(i);
			double aObs = data.getA(i);
			

			
			// Expected v under model
			double vExpected = this.getExpectedV(aObs);
			
			
			// Error
			if (logInput.get()) {
				error = Math.log(vExpected) - Math.log(vObs);
				
				// Cannot log zero, but the expected value should be zero anyway
				if (vExpected == 0) {
					if (vObs != 0) {
						logP = Double.NEGATIVE_INFINITY;
						return logP;
					}else {
						continue;
					}
				}
				
				
			}else {
				error = vExpected - vObs;
			}
			
			
			// Probability density of error under normal distribution
			logP +=  dist.logDensity(error);
			 
			
		}
		
		
		return logP;
	}
	
	
	/**
	 * Prepare for a series of 'getExpectedV' calculations from the current state
	 * Return false if state is illegal
	 */
	public boolean prepare() {
		return true;
	}
	
	
	/**
	 * Get the active Hill coefficient
	 * @return
	 */
	public double getHill() {
		
		if (hillIndicatorInput.get() == null) return 1;
		
		
		// Regular
		if (hillIndicatorInput.get().getValue() == -1) {
			return hillInput.get().getValue();
		}
		
		// Inversion
		if (hillIndicatorInput.get().getValue() == 1) {
			return 1.0 / hillInput.get().getValue();
		}
		
		return 1;
	}
	
	
	
	/**
	 * The expected rate velocity under the model, as a function of reactant concentration
	 * @param a
	 * @return
	 */
	public double getExpectedV(double a) {
		double km = kmInput.get().getValue();
		double vmax = vmaxInput.get().getValue();
		double hill = getHill();
		if (hill != 1) {
			a = Math.pow(a, hill);
			km = Math.pow(km, hill);
		}
		
		return (vmax*a) / (km + a);
	}
	
	
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}
	
	

	@Override
	public void init(PrintStream out) {
		KineticData data = dataInput.get();
		out.print("hill\t");
		for (int i = 0; i < data.getNumberOfObservations(); i ++) {
			out.print("a." + (i+1)  + "\t");
			if (data.hasV(i)) out.print("v.obs." + (i+1)  + "\t");
			out.print("v.pred." + (i+1)  + "\t");
			out.print("v.sample." + (i+1)  + "\t");
		}
		
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		KineticData data = dataInput.get();
		this.prepare();
		out.print(this.getHill() + "\t");
		for (int i = 0; i < data.getNumberOfObservations(); i ++) {
			
			// Sample a
			double a = data.getObservation(i).getA();
			
			// Calculate v 
			double v = this.getExpectedV(a);
			
			
			// Sample v with error
			double sd = sdInput.get().getValue();
			double vErr;
			if (logInput.get()) {
				vErr = Math.exp(Math.log(v) + Randomizer.nextGaussian() * sd);
			}else {
				vErr = v + Randomizer.nextGaussian() * sd;
			}
			
			out.print(a  + "\t");
			if (data.hasV(i)) out.print(data.getObservation(i).getV()  + "\t");
			out.print(v  + "\t");
			out.print(vErr  + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}

}

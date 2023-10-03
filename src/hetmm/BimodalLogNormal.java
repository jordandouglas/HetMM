package hetmm;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math3.distribution.NormalDistribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;


@Description("A bimodal log-normal distribution")
public class BimodalLogNormal extends Distribution {
	
	public final Input<RealParameter> xInput = new Input<>("x", "parameter", Input.Validate.REQUIRED);
	public final Input<RealParameter> m1Input = new Input<>("m1", "mean 1", Input.Validate.REQUIRED);
	public final Input<RealParameter> m2Input = new Input<>("m2", "mean 2", Input.Validate.REQUIRED);
	public final Input<RealParameter> s1Input = new Input<>("s1", "std dev 1", Input.Validate.REQUIRED);
	public final Input<RealParameter> s2Input = new Input<>("s2", "std dev 2", Input.Validate.REQUIRED);
	
	
	@Override
	public double calculateLogP() {
		logP = 0;
		
		double m1 = m1Input.get().getValue();
		double m2 = m2Input.get().getValue();
		double s1 = s1Input.get().getValue();
		double s2 = s2Input.get().getValue();
		NormalDistribution d1 = new NormalDistribution(m1, s1);
		NormalDistribution d2 = new NormalDistribution(m2, s2);
		
		
		
		for (int i = 0; i < xInput.get().getDimension(); i++) {
			double x = xInput.get().getValue(i);
			if (x <= 0) {
				logP = Double.NEGATIVE_INFINITY;
				return logP;
			}
			
			double logx = Math.log(x);
			//double p1 = d1.logDensity(logx);
			//double p2 = d2.logDensity(logx);
			
			double p1 = 0.5 * d1.density(logx);
			double p2 = 0.5 * d2.density(logx);
			
			logP += Math.log(p1) + Math.log(p2);
			//Log.warning(p1 + " " + p2 + " " + x + " " + logx);
		}
		
		
		
		return logP;
	}
	



	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(xInput.get().getID());
		return args;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(m1Input.get().getID());
		conds.add(m2Input.get().getID());
		conds.add(s1Input.get().getID());
		conds.add(s2Input.get().getID());
		return conds;
	}

	@Override
	public void sample(State state, Random random) {
		
		if (this.sampledFlag) return;
		sampleConditions(state, random);
		this.sampledFlag = true;
		
		double m1 = m1Input.get().getValue();
		double m2 = m2Input.get().getValue();
		double s1 = s1Input.get().getValue();
		double s2 = s2Input.get().getValue();
		NormalDistribution d1 = new NormalDistribution(m1, s1);
		NormalDistribution d2 = new NormalDistribution(m2, s2);
		
		for (int i = 0; i < xInput.get().getDimension(); i++) {
			
			// Sample a distribution
			double x = 0;
			try {
			boolean useD1 = random.nextBoolean();
			if (useD1) {
				x = d1.inverseCumulativeProbability(random.nextFloat());
			}else {
				x = d2.inverseCumulativeProbability(random.nextFloat());
			}
			} catch(Exception e) {
				
			}
			
			x = Math.exp(x);
			xInput.get().setValue(i, x);
			
		}
		
		
	}

	
	
}

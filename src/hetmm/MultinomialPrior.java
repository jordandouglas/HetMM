package hetmm;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("A prior on an integer - 1 dimensional only")
public class MultinomialPrior extends Distribution {
	
	
	final public Input<RealParameter> pInput = new Input<>("p", "probabiltiy of each element being that index", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> xInput = new Input<>("x", "", Input.Validate.REQUIRED);
	
	
	@Override
	public void initAndValidate() {
		
		if (xInput.get().getUpper() - xInput.get().getLower() + 1 != pInput.get().getDimension()) {
			throw new IllegalArgumentException("Mismatching dimension between p and the range of x");
		}
	}
	
	public double calculateLogP() {
		
        logP = 0;
        
        IntegerParameter x = xInput.get();

        int val = x.getValue(0);
        for (int i = 0; i < pInput.get().getDimension(); i++) {
        	
        	int valOfIndex = xInput.get().getLower() + i;
        	if (val == valOfIndex) {
        		double pVal = pInput.get().getValue(i);
        		logP = Math.log(pVal);
        		return logP;
        	}
        	
        	
        }
        
        
        logP = Double.NEGATIVE_INFINITY;
        return logP;
	}
	
	@Override
	public List<String> getArguments() {
		List<String> conds = new ArrayList<>();
		conds.add(xInput.get().getID());
		return conds;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(pInput.get().getID());
		return conds;
	}

	@Override
	public void sample(State state, Random random) {
		
		
		if (this.sampledFlag) return;
		this.sampleConditions(state, random);
		IntegerParameter x = xInput.get();
        double[] p = pInput.get().getDoubleValues();
        int n = x.getDimension();
        
        int val = Randomizer.randomChoicePDF(p) + xInput.get().getLower();
        x.setValue(0, val);
		this.sampledFlag = true;
		
	}

}

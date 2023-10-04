package hetmm;


import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;

@Description("Inverts a number")
public class Inversion extends CalculationNode implements Function {
	
	
	final public Input<Function> xInput = new Input<>("x", "function to invert", Input.Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}


	@Override
	public int getDimension() {
		return xInput.get().getDimension();
	}

	@Override
	public double getArrayValue(int dim) {
		double x = xInput.get().getArrayValue(dim);
		return 1.0 / x;
	}

}

namespace ISAAR.MSolve.Optimization.Problems
{
    /// <summary>
    /// Represents a specific design, i.e. a configuration of the design variables. In order to input the objective and 
    /// contraint functions of the optimization problem, client code must: a) implement this interface, b) evaluate 
    /// those functions in private, c) assign the calculated values in the corresponding properties.
    /// </summary>
    /// <remarks>This interface simplifies the evaluation of various functions (objective or constraint) that have the 
    /// same arguments (design variables) and require common steps, i.e. a FEM simulation. Thus the common steps are 
    /// executed only once in the constructor of the concrete classes implementing <see cref="IDesign"/>. Client code 
    /// is also repsonsible for evaluating these functions in the constructor and assigning them to the corresponding 
    /// properties. Finally the optimization framework can access the calculated values as required.</remarks>
    public interface IDesign
    {
        double[] ObjectiveValues { get; }
        double[] ConstraintValues { get; }
    }

    /// <summary>
    /// A Factory interface to create instances of <see cref="IDesign"/>. Client code must implement this interface to 
    /// provide the optimization framework with a way to instantiate <see cref="IDesign"/> objects.
    /// </summary>
    public interface IDesignFactory
    {
        /// <summary>
        /// Creates an object implementing <see cref="IDesign"/>
        /// </summary>
        /// <param name="x">The design variables specifying the requested design</param>
        /// <returns></returns>
        IDesign CreateDesign(double[] x);
    }
}
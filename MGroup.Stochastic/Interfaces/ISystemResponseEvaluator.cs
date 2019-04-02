namespace MGroup.Stochastic.Interfaces
{
    /// <summary>The interface of the classes that evaluate the selected response for each iteration of the stochastic process</summary>
    public interface ISystemResponseEvaluator
    {
        double[] Evaluate(int iteration);
    }
}

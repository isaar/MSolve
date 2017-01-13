namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IStochasticCoefficientsProvider
    {
        double GetCoefficient(double meanValue, double[] coordinates);
    }
}

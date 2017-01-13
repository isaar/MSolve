namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IStochasticMaterialCoefficientsProvider
    {
        double[] RandomVariables { get; set; }
        double GetCoefficient(double meanValue, double[] coordinates);
    }
}

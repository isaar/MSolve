namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IStochasticFiniteElement : IFiniteElement
    {
        IStochasticCoefficientsProvider CoefficientsProvider { get; set; }
    }
}

namespace MGroup.Stochastic.Interfaces
{
    /// <summary>The interface of the classes that generate random realizations of the stochastic system parameter</summary>
    public interface IUncertainParameterRealizer
    {
        double Realize(int iteration, IStochasticDomainMapper domainMapper, double[] parameters);
    }
}

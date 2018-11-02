namespace MGroup.Stochastic.Interfaces
{
    public interface IUncertainParameterRealizer
    {
        double Realize(int iteration, IStochasticDomainMapper domainMapper, double[] parameters);
    }
}

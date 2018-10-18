namespace MGroup.Stochastic.Interfaces
{
    public interface IUncertainParameterRealizer
    {
        double[] Realize(int iteration, double[] parameters);
    }
}

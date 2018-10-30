namespace MGroup.Stochastic.Interfaces
{
    public interface IUncertainParameterRealizer
    {
        double Realize(int iteration, int parameters);
    }
}

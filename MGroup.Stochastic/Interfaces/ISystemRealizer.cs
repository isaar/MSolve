namespace MGroup.Stochastic.Interfaces
{
    /// <summary>The interface of the classes that create new model realizations</summary>
    public interface ISystemRealizer
    {
        void Realize(int iteration);
    }
}


//TODO: this class should be used instead of INode in solvers, analyzers, etc. Since those classes do not care if the problem is 
//      2D or 3D, the coordinates will not be used.
namespace ISAAR.MSolve.Discretization
{
    /// <summary>
    /// Tagging interface for discrete point where interesting quantities are known or computed. E.g. in FEM such points are 
    /// called nodes.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IDiscretePoint
    {
        int ID { get; }
    }
}

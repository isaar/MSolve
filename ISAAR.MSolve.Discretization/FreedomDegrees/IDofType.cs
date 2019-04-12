//WARNING: DofTable depends on this being immutable.
//TODO: The implementing classes (structural, thermal, pore pressure) should be moved to projects that correspond to that 
//      differential equation.  
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// Tagging interface for degrees of freedom. Concrete implementations must be immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface IDofType
    {
    }
}

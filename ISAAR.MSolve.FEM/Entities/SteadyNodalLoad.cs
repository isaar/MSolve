using ISAAR.MSolve.Discretization.FreedomDegrees;

//TODO: This is probably covered by Load.cs
namespace ISAAR.MSolve.FEM.Entities
{
    /// <summary>
    /// You should use <see cref="Load"/> instead.
    /// </summary>
    public class SteadyNodalLoad: ITimeDependentNodalLoad
    {
        private readonly double constantloadAmount;

        public SteadyNodalLoad(double constantloadAmount)
        {
            this.constantloadAmount = constantloadAmount;
        }

        public Node Node { get; set; }
        public IDofType DOF { get; set; }

        public double GetLoadAmount(int timeStep)
        {
            return constantloadAmount;
        }
    }
}

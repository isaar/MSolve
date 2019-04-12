using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.IGA.Entities.Loads
{
    public class ControlPointLoad : Load
    {
        public ControlPoint ControlPoint { get; set; }
        public IDofType DOF { get; set; }
        public double Amount { get; set; }
    }
}

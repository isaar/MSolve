using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Discretization
{
    public class Constraint
    {
        //public Node Node { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}


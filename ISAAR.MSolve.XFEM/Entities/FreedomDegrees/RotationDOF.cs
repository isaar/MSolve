using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    public class RotationDOF: IDOF
    {
        public static readonly RotationDOF X = new RotationDOF("Rotation around X axis");
        public static readonly RotationDOF Y = new RotationDOF("Rotation around Y axis");
        public static readonly RotationDOF Z = new RotationDOF("Rotation around Z axis");

        private readonly string name;

        private RotationDOF(string name) // No more displacement dofs can be created
        {
            this.name = name;
        }

        public override string ToString()
        {
            return name;
        }
    }
}

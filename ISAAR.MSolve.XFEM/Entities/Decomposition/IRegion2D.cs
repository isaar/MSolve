using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.Elements;

// TODO: Should this be here or in Geometry?
namespace ISAAR.MSolve.XFEM.Entities.Decomposition
{
    public enum NodePosition
    {
        Internal, Boundary, External
    }

    interface IRegion2D
    {
        NodePosition FindRelativePosition(XNode2D node);
    }
}

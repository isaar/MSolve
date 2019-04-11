using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ICollocationElement:IElement_v2
    {
        INode CollocationPoint { get; set; }
    }
}

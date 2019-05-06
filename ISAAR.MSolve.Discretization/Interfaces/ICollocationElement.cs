using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ICollocationElement:IElement
    {
        INode CollocationPoint { get; set; }
    }
}

using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface INode
    {
		int ID { get; set; }
		double X { get; set; }
		double Y { get; set; }
		double Z { get; set; }
	}
}

using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
	/// <summary>
	/// Vertex of a finite element in a 3-dimensional space. Immutable.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class Node3D : Node, IComparable<Node3D>
	{
		/// <summary>
		/// Instantiates a <see cref="Node3D"/>
		/// </summary>
		/// <param name="id">A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node2D"/> 
		///     from the rest. <see cref="ID"/> &gt;= 0.</param>
		/// <param name="x">The coordinate of the point along axis X.</param>
		/// <param name="y">The coordinate of the point along axis Y.</param>
		/// <param name="z">The coordinate of the point along axis Z.</param>
		public Node3D(int id, double x, double y, double z)
		{
			if (id<0) throw new ArgumentException("The parameter id must be non negative, but was: "+id);
			this.ID = id;
			this.X = x;
			this.Y = y;
			this.Z = z;
		}

		public int CompareTo(Node3D other)
		{
			return this.ID - other.ID;
		}

	}
}

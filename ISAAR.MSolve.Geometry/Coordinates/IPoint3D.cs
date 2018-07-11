namespace ISAAR.MSolve.Geometry.Coordinates
{
	/// <summary>
	/// Point in a 3-dimensional space. Immutable
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public interface IPoint3D
	{
		/// <summary>
		/// Vector with the coordinates of the point. Length = 3.
		/// </summary>
		double[] Coordinates { get; }
	}
}
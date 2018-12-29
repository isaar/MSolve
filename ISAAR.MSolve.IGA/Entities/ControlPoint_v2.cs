using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;

namespace ISAAR.MSolve.IGA.Entities
{
	public class ControlPoint_v2 : ControlPoint
	{
		private readonly List<Constraint> constraints = new List<Constraint>();
		private readonly Dictionary<int, Patch_v2> patchesDictionary_v2 = new Dictionary<int, Patch_v2>();

		public List<Constraint> Constraints_v2 => constraints;

		public Dictionary<int, Patch_v2> PatchesDictionary_v2 => patchesDictionary_v2;

		public void BuildPatchesDictionary_v2()
		{
			foreach (Element element in elementsDictionary.Values)
				if (!patchesDictionary_v2.ContainsKey(element.Patch_v2.ID))
					patchesDictionary_v2.Add(element.Patch_v2.ID, element.Patch_v2);
		}
	}
}
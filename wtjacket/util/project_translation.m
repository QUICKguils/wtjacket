function qProjected = project_translation(q, direction, nodeList, nodeLabel)
	% PROJECT_TRANSLATION  Translation from a displacement.
	%
	% This function extracts the translation motion of a displacement at
	% a given node label, and projects it in a given direction.
	%
	% Arguments:
	%	q         (nDofFreexnTime double) -- Displacements [m].
	%	direction (1x3 double)            -- Direction of projection [m].
	%	nodeList  {1xN Node}              -- Cell list of nodes.
	%	nodeLabel (int)                   -- Label of the node to study.
	% Return:
	%	qProjected (1xnTime double) -- Projected translation [m].

	qX = q(nodeList{nodeLabel}.dof(1), :) * direction(1);
	qY = q(nodeList{nodeLabel}.dof(2), :) * direction(2);
	qZ = q(nodeList{nodeLabel}.dof(3), :) * direction(3);

	qProjected = qX + qY + qZ;
end
struct PolygonalObstacle
    center::Vector
    dx::Real
    dy::Real
    dz::Real
    normals::Array{Real,2}
    corners::Array{Real,2}
end
function PolygonalObstacle(center::Vector,width::Vector)
  dx = width[1]/2.0
  dy = width[2]/2.0
  dz = width[3]/2.0
  corners = set_corners(center, dx, dy, dz)
  normals = set_normals(center, dx, dy, dz)
  return PolygonalObstacle(center,dx,dy,dz,normals,corners)
end

function ISS_koz(kozN::Integer)
  btms_lft, tops_rgt = repeat([zeros(3)], kozN), repeat([zeros(3)], kozN)
  	btms_lft[1],  tops_rgt[1]  = Vector([ 11.9,2.7,0.0 ]), Vector([ 12.0,12.0,6.0 ])  		# 1
	btms_lft[2],  tops_rgt[2]  = Vector([ 11.6,1.2,0.0 ]), Vector([ 12.0,2.7,6.0 ])  		# 2
	btms_lft[3],  tops_rgt[3]  = Vector([ 11.6,0.8,0.0 ]), Vector([ 12.0,1.2,6.0 ])  		# 3
	btms_lft[4],  tops_rgt[4]  = Vector([ 11.6,-12.0,0.0 ]), Vector([ 12.0,-0.8,6.0 ])  	# 4
	btms_lft[5],  tops_rgt[5]  = Vector([ -12.0,-12.0,0.0 ]), Vector([ 11.6,-2.7,6.0 ])  	# 6
	btms_lft[6],  tops_rgt[6]  = Vector([ -12.0,-2.7,0.0 ]), Vector([ 10.3,-1.2,6.0 ])  	# 7
	btms_lft[7],  tops_rgt[7]  = Vector([ -12.0,-1.2,0.0 ]), Vector([ 7.7,-0.6,6.0 ])  		# 8
	btms_lft[8],  tops_rgt[8]  = Vector([ -12.0,-0.6,0.0 ]), Vector([ 5.9,12.0,6.0 ])  		# 9
	btms_lft[9],  tops_rgt[9]  = Vector([ 5.9,0.6,0.0 ]), Vector([ 7.7,12.0,6.0 ])  		# 10
	btms_lft[10],  tops_rgt[10]  = Vector([ 9.6,7.3,0.0 ]), Vector([ 11.9,12.0,6.0 ])  		# 11
	btms_lft[11],  tops_rgt[11]  = Vector([ 7.7,2.7,0.0 ]), Vector([ 9.6,12.0,6.0 ])  		# 12
	btms_lft[12],  tops_rgt[12]  = Vector([ 7.7,1.2,0.0 ]), Vector([ 10.2,2.7,6.0 ])  		# 13
	#
	# btms_lft[14],  tops_rgt[14]  = Vector([ 11.9,7.3,0.0 ]), Vector([ 9.6,2.7,3.8 ])  		# 1b
	# btms_lft[15],  tops_rgt[15]  = Vector([ 11.9,7.3,5.9 ]), Vector([ 9.6,2.7,6.0 ])  		# 1t
	# btms_lft[16],  tops_rgt[16]  = Vector([ 11.6,2.7,0.0 ]), Vector([ 10.2,1.2,4.2 ])		# 2b
	# btms_lft[17],  tops_rgt[17]  = Vector([ 11.6,2.7,5.5 ]), Vector([ 10.2,1.2,6.0 ])		# 2t
	# btms_lft[18],  tops_rgt[18]  = Vector([ 12.0,0.8,0.0 ]), Vector([ 11.6,-0.8,4.1 ])		# 3b
	# btms_lft[19],  tops_rgt[19]  = Vector([ 12.0,0.8,5.5 ]), Vector([ 11.6,-0.8,6.0 ])		# 3t
	# btms_lft[20],  tops_rgt[20]  = Vector([ 11.6,-1.2,0.0 ]), Vector([ 10.3,-2.7,4.3 ])		# 4b
	# btms_lft[21],  tops_rgt[21]  = Vector([ 11.6,-1.2,5.4 ]), Vector([ 10.3,-2.7,6.0 ])		# 4t
	# btms_lft[22],  tops_rgt[22]  = Vector([ 11.6,1.2,0.0 ]), Vector([ 7.7,-1.2,3.7 ])		# 5b
	# btms_lft[23],  tops_rgt[23]  = Vector([ 11.6,1.2,6.0 ]), Vector([ 7.7,-1.2,6.0 ])		# 5t
	# btms_lft[24],  tops_rgt[24]  = Vector([ 7.7,0.6,0.0 ]), Vector([ 5.9,-0.6,4.3 ])		# 6b
	# btms_lft[25],  tops_rgt[25]  = Vector([ 7.7,0.6,5.4 ]), Vector([ 5.9,-0.6,6.0 ])		# 6t


	# btms_lft[1],  tops_rgt[1]  = Vector([ 3.0,-0.6, 4.2]), Vector([ 5.9, 0.6, 5.4])  # 1
	# btms_lft[2],  tops_rgt[2]  = Vector([ 5.0, 0.6, 3.7]), Vector([ 7.7, 1.2, 6.0])  # 2
	# btms_lft[3],  tops_rgt[3]  = Vector([ 5.9, 1.2, 3.7]), Vector([10.2, 2.7, 6.0])  # 3
	# btms_lft[4],  tops_rgt[4]  = Vector([ 8.0, 2.7, 3.8]), Vector([ 9.6, 7.3, 5.9])  # 4
	# btms_lft[5],  tops_rgt[5]  = Vector([ 9.6, 7.3, 3.8]), Vector([11.9, 9.0, 5.9])  # 5
	# btms_lft[6],  tops_rgt[6]  = Vector([11.9, 2.7, 3.8]), Vector([13.0, 7.3, 5.9])  # 6
	# btms_lft[7],  tops_rgt[7]  = Vector([11.6, 1.2, 3.8]), Vector([12.0, 2.7, 5.9])  # 7
	# btms_lft[8],  tops_rgt[8]  = Vector([11.6, 0.8, 3.7]), Vector([12.0, 1.2, 6.0])  # 8
	# btms_lft[9],  tops_rgt[9]  = Vector([12.0,-0.8, 3.7]), Vector([13.0, 0.8, 6.0])  # 9
	# btms_lft[10], tops_rgt[10] = Vector([11.6,-1.2, 3.7]), Vector([12.0,-0.8, 6.0])  # 10
	# btms_lft[11], tops_rgt[11] = Vector([11.6,-2.7, 4.3]), Vector([12.0,-1.2, 5.4])  # 11
	# btms_lft[12], tops_rgt[12] = Vector([10.3,-4.0, 4.3]), Vector([11.6,-2.7, 5.4])  # 12
	# btms_lft[13], tops_rgt[13] = Vector([ 5.9,-4.0, 0.0]), Vector([10.3,-1.2, 0.0])  # 13
	# btms_lft[14], tops_rgt[14] = Vector([ 5.9,-1.2, 3.7]), Vector([ 7.7,-0.6, 6.0])  # 14
	# btms_lft[15], tops_rgt[15] = Vector([ 5.9,-1.2, 5.4]), Vector([ 7.7, 1.2, 6.0])  # 15
	# btms_lft[16], tops_rgt[16] = Vector([ 5.9,-1.2, 3.0]), Vector([ 7.7, 1.2, 4.2])  # 16
	# btms_lft[17], tops_rgt[17] = Vector([10.2, 1.2, 5.5]), Vector([11.6, 2.7, 7.0])  # 17
	# btms_lft[18], tops_rgt[18] = Vector([10.2, 1.2, 3.0]), Vector([11.6, 2.7, 4.2])  # 18
	# btms_lft[19], tops_rgt[19] = Vector([ 9.6, 2.7, 3.0]), Vector([11.9, 7.3, 3.8])  # 19
	# btms_lft[20], tops_rgt[20] = Vector([ 9.6, 2.7, 5.9]), Vector([11.9, 7.3, 7.0])  # 20
	# btms_lft[21], tops_rgt[21] = Vector([10.3,-2.7, 5.4]), Vector([11.6,-1.2, 7.0])  # 21
	# btms_lft[22], tops_rgt[22] = Vector([10.3,-2.7, 3.0]), Vector([11.6,-1.2, 4.3])  # 22
	# btms_lft[23], tops_rgt[23] = Vector([ 7.7,-1.2, 3.0]), Vector([11.6, 1.2, 3.7])  # 23
	# btms_lft[24], tops_rgt[24] = Vector([ 7.7,-1.2, 6.0]), Vector([11.6, 1.2, 7.0])  # 24
	# btms_lft[25], tops_rgt[25] = Vector([11.6,-0.8, 3.0]), Vector([12.0, 0.8, 4.1])  # 25
	# btms_lft[26], tops_rgt[26] = Vector([11.6,-0.8, 5.5]), Vector([12.0, 0.8, 7.0])  # 26
	keepout_zones = []
	for (btm, top) in zip(btms_lft, tops_rgt)
	    center, width = (top+btm)/2., (top-btm)
	    append!(keepout_zones, [PolygonalObstacle(center,width)] )
	end
	return keepout_zones
end

function set_normals(center::Vector,dx::Real,dy::Real,dz::Real)
  # faces are numbered using:
  #  1 is the +z face (top)
  #  2 is the +x face (front)
  #  3 is the +y face (right)
  #  4 is the -x face (rear)
  #  5 is the -y face (left)
  #  6 is the -z face (bottom)
  n1 = Vector([ 0.,  0.,  1.])
  n2 = Vector([ 1.,  0.,  0.])
  n3 = Vector([ 0.,  1.,  0.])
  n4 = Vector([-1.,  0.,  0.])
  n5 = Vector([ 0., -1.,  0.])
  n6 = Vector([ 0.,  0., -1.])
  return hcat(n1, n2, n3, n4, n5, n6)
end
function get_normals(obs::PolygonalObstacle)
  n = obs.normals
  return n[:,1], n[:,2], n[:,3], n[:,4], n[:,5], n[:,6]
end

function set_corners(center::Vector,dx::Real, dy::Real, dz::Real)
  # corners are numbered using the three faces that meet there.
  # faces are numbered using:
  #  1 is the +z face (top)
  #  2 is the +x face (front)
  #  3 is the +y face (right)
  #  4 is the -x face (rear)
  #  5 is the -y face (left)
  #  6 is the -z face (bottom)
  cx, cy, cz = center[1], center[2], center[3]
  c125 = Vector([cx + dx, cy + dy, cz + dz])
  c256 = Vector([cx + dx, cy + dy, cz - dz])
  c236 = Vector([cx + dx, cy - dy, cz - dz])
  c123 = Vector([cx + dx, cy - dy, cz + dz])
  c134 = Vector([cx - dx, cy - dy, cz + dz])
  c346 = Vector([cx - dx, cy - dy, cz - dz])
  c456 = Vector([cx - dx, cy + dy, cz - dz])
  c145 = Vector([cx - dx, cy + dy, cz + dz])
  return hcat(c125, c256, c236, c123, c134, c346, c456, c145)
end
function get_corners(obs::PolygonalObstacle)
  c = obs.corners
  c125, c256, c236, c123, c134, c346, c456, c145 = (c[:,1],
        c[:,2], c[:,3], c[:,4], c[:,5], c[:,6], c[:,7], c[:,8])
  return c125, c256, c236, c123, c134, c346, c456, c145
end

function is_in_halfspace(r::Vector,i::Integer,obs::PolygonalObstacle)
  # check to see if the vector r lies in the halfspace
  # defined by the i-th outward facing normal of the obstacle obs.
  # faces are numbered using:
  #  1 is the +z face (top)
  #  2 is the +x face (front)
  #  3 is the -y face (left)
  #  4 is the -x face (rear)
  #  5 is the +y face (right)
  #  6 is the -z face (bottom)
  if i>6
    error("Input i must be no greater than 6")
  end
  c = obs.center
  n = obs.normals[:,i]

  if (i == 2 || i == 4)
    return (dot(r-c,n) >= obs.dx)
  elseif (i == 3 || i == 5)
    return (dot(r-c,n) >= obs.dy)
  elseif (i == 1 || i == 6)
    return (dot(r-c,n) >= obs.dz)
  end
end

function signed_distance(r::Vector,obs::PolygonalObstacle)
  cntr = obs.center
  cx, cy, cz = cntr[1], cntr[2], cntr[3]
  dx, dy, dz = obs.dx, obs.dy, obs.dz
  c125, c256, c236, c123, c134, c346, c456, c145 = get_corners(obs)
  n1, n2, n3, n4, n5, n6 = get_normals(obs)

  # check in which halfspace(s) vector r is in.
  # hlfsps[k] is true if x is in k-th halfspace
  hlfsps = falses(6)
  for k = 1:6
    hlfsps[k] = is_in_halfspace(r,k,obs)
  end

  # check if r is inside the obstacle
  inside = !hlfsps[1] && !hlfsps[2] && !hlfsps[3] && !hlfsps[4] && !hlfsps[5] && !hlfsps[6]

  if inside
    # inward facing normals are -ve of the outward facing normals
    projs = zeros(6)
    projs[1] = dot(r-cntr,-n1)
    projs[2] = dot(r-cntr,-n2)
    projs[3] = dot(r-cntr,-n3)
    projs[4] = dot(r-cntr,-n4)
    projs[5] = dot(r-cntr,-n5)
    projs[6] = dot(r-cntr,-n6)
    id = argmin(projs)
    if id==1
      close_pt = [ r[1]; r[2]; cz+dz ]
    elseif id==2
      close_pt = [ cx+dx; r[2]; r[3] ]
    elseif id==3
      close_pt = [ r[1]; cy-dy; r[3] ]
    elseif id==4
      close_pt = [ cx-dx; r[2]; r[3] ]
    elseif id==5
      close_pt = [ r[1]; cy+dy; r[3] ]
    else
      close_pt = [ r[1]; r[2]; cz-dz ]
    end
    # since we're inside the obstacle, signed distance returns
    # a negative number
    return -norm(r-close_pt), -close_pt
  else
    # r is outside a corner
    if (hlfsps[1] && hlfsps[2] && hlfsps[5])
      close_pt = c125
    elseif (hlfsps[2] && hlfsps[5] && hlfsps[6])
      close_pt = c256
    elseif (hlfsps[2] && hlfsps[3] && hlfsps[6])
      close_pt = c236
    elseif (hlfsps[1] && hlfsps[2] && hlfsps[3])
      close_pt = c123
    elseif (hlfsps[1] && hlfsps[3] && hlfsps[4])
      close_pt = c134
    elseif (hlfsps[3] && hlfsps[4] && hlfsps[6])
      close_pt = c346
    elseif (hlfsps[4] && hlfsps[5] && hlfsps[6])
      close_pt = c456
    elseif (hlfsps[1] && hlfsps[4] && hlfsps[5])
      close_pt = c145
    # r is outside an edge
    elseif (hlfsps[1] && hlfsps[2])
      close_pt = [ cx+dx; r[2]; cz+dz ]
    elseif (hlfsps[1] && hlfsps[5])
      close_pt = [ r[1]; cy+dy; cz+dz ]
    elseif (hlfsps[1] && hlfsps[4])
      close_pt = [ cx-dx; r[2]; cz+dz ]
    elseif (hlfsps[1] && hlfsps[3])
      close_pt = [ r[1]; cy-dy; cz+dz ]
    elseif (hlfsps[2] && hlfsps[5])
      close_pt = [ cx+dx; cy+dy; r[3] ]
    elseif (hlfsps[2] && hlfsps[3])
      close_pt = [ cx+dx; cy-dy; r[3] ]
    elseif (hlfsps[3] && hlfsps[4])
      close_pt = [ cx-dx; cy-dy; r[3] ]
    elseif (hlfsps[4] && hlfsps[5])
      close_pt = [ cx-dx; cy+dy; r[3] ]
    elseif (hlfsps[2] && hlfsps[6])
      close_pt = [ cx+dx; r[2]; cz-dz ]
    elseif (hlfsps[3] && hlfsps[6])
      close_pt = [ r[1]; cy-dy; cz-dz ]
    elseif (hlfsps[4] && hlfsps[6])
      close_pt = [ cx-dx; r[2]; cz-dz ]
    elseif (hlfsps[5] && hlfsps[6])
      close_pt = [ r[1]; cy+dy; cz-dz ]
    # r is outside one face
    elseif hlfsps[1]
      close_pt = [ r[1]; r[2]; cz+dz ]
    elseif hlfsps[2]
      close_pt = [ cx+dx; r[2]; r[3] ]
    elseif hlfsps[3]
      close_pt = [ r[1]; cy-dy; r[3] ]
    elseif hlfsps[4]
      close_pt = [ cx-dx; r[2]; r[3] ]
    elseif hlfsps[5]
      close_pt = [ r[1]; cy+dy; r[3] ]
    elseif hlfsps[6]
      close_pt = [ r[1]; r[2]; cz-dz ]
    end
    # since we're outside the obstacle, signed distance returns
    # a positive number (first output)
    return norm(r-close_pt), close_pt
  end
end

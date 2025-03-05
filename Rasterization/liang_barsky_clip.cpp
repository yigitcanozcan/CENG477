if(currentMesh->type==0){
			OurVertex v0 = currentTriangle.vertices[0];
            OurVertex v1 = currentTriangle.vertices[1];
            OurVertex v2 = currentTriangle.vertices[2];
            Color c0 = Color (currentTriangle.vertices[0].color.x, currentTriangle.vertices[0].color.y, currentTriangle.vertices[0].color.z);
            Color c1 = Color (currentTriangle.vertices[1].color.x, currentTriangle.vertices[1].color.y, currentTriangle.vertices[1].color.z);
            Color c2 = Color (currentTriangle.vertices[2].color.x, currentTriangle.vertices[2].color.y, currentTriangle.vertices[2].color.z);


                std::pair <Vec4, Color> L01_pair1 = std::make_pair(position0, c0);
                std::pair <Vec4, Color> L01_pair2 = std::make_pair(position1, c1);

                // Line-2 => L12
                std::pair <Vec4, Color> L12_pair1 = std::make_pair(position1, c1);
                std::pair <Vec4, Color> L12_pair2 = std::make_pair(position2, c2);

                // Line-3 => L20
                std::pair <Vec4, Color> L20_pair1 = std::make_pair(position2, c2);
                std::pair <Vec4, Color> L20_pair2 = std::make_pair(position0, c0);

                /* Clipping Phase */
                bool L01_visibility = clipLine(L01_pair1, L01_pair2);
                bool L12_visibility = clipLine(L12_pair1, L12_pair2);
                bool L20_visibility = clipLine(L20_pair1, L20_pair2);

                /* Viewport Transformation Phase */
                // L01
                L01_pair1.first = multiplyMatrixWithVec4(viewportTransformationMatrix, L01_pair1.first);
                L01_pair2.first = multiplyMatrixWithVec4(viewportTransformationMatrix, L01_pair2.first);

                // L12
                L12_pair1.first = multiplyMatrixWithVec4(viewportTransformationMatrix, L12_pair1.first);
                L12_pair2.first = multiplyMatrixWithVec4(viewportTransformationMatrix, L12_pair2.first);

                // L20
                L20_pair1.first = multiplyMatrixWithVec4(viewportTransformationMatrix, L20_pair1.first);
                L20_pair2.first = multiplyMatrixWithVec4(viewportTransformationMatrix, L20_pair2.first);

				OurVertex L01_first;
				L01_first.position.x = L01_pair1.first.x;
				L01_first.position.y = L01_pair1.first.y;
				L01_first.position.z = L01_pair1.first.z;
				L01_first.color.x = L01_pair1.second.r;
				L01_first.color.y = L01_pair1.second.g;
				L01_first.color.z = L01_pair1.second.b;
				OurVertex L01_second;
				L01_second.position.x = L01_pair2.first.x;
				L01_second.position.y = L01_pair2.first.y;
				L01_second.position.z = L01_pair2.first.z;
				L01_second.color.x = L01_pair2.second.r;
				L01_second.color.y = L01_pair2.second.g;
				L01_second.color.z = L01_pair2.second.b;

				OurVertex L12_first;
				L12_first.position.x = L12_pair1.first.x;
				L12_first.position.y = L12_pair1.first.y;
				L12_first.position.z = L12_pair1.first.z;
				L12_first.color.x = L12_pair1.second.r;
				L12_first.color.y = L12_pair1.second.g;
				L12_first.color.z = L12_pair1.second.b;
				OurVertex L12_second;
				L12_second.position.x = L12_pair2.first.x;
				L12_second.position.y = L12_pair2.first.y;
				L12_second.position.z = L12_pair2.first.z;
				L12_second.color.x = L12_pair2.second.r;
				L12_second.color.y = L12_pair2.second.g;
				L12_second.color.z = L12_pair2.second.b;

				OurVertex L20_first;
				L20_first.position.x = L20_pair1.first.x;
				L20_first.position.y = L20_pair1.first.y;
				L20_first.position.z = L20_pair1.first.z;
				L20_first.color.x = L20_pair1.second.r;
				L20_first.color.y = L20_pair1.second.g;
				L20_first.color.z = L20_pair1.second.b;
				OurVertex L20_second;
				L20_second.position.x = L20_pair2.first.x;
				L20_second.position.y = L20_pair2.first.y;
				L20_second.position.z = L20_pair2.first.z;
				L20_second.color.x = L20_pair2.second.r;
				L20_second.color.y = L20_pair2.second.g;
				L20_second.color.z = L20_pair2.second.b;

                /* Final step of FRP - Rasterize the Line and fill the image for this model */
                if(L01_visibility)
                    drawLine(L01_first, L01_second, image); // L01
                if(L12_visibility)
                    drawLine(L12_first, L12_second, image); // L12
                if(L20_visibility)
                    drawLine(L20_first, L20_second, image); // L20
			}
			
			
			

//liang-barsky
bool visible(double den, double num, double & t_E, double & t_L) {
    double t = num / den;
    if (den > 0) {
        //PE
        if (t > t_L)
            return false;
        else if (t > t_E)
            t_E = t;
    }
    else if (den < 0) {
        //PL
        if (t < t_E)
            return false;
        else if (t < t_L)
            t_L = t;
    }
    else if (num > 0) {
        //parallel
        return false;
    }
    return true;
}

bool clipLine(std::pair<Vec4, Color> & pair1, std::pair<Vec4, Color> & pair2) {
    bool isVisible = false;
    Vec4 v0 = pair1.first, v1 = pair2.first;
    Color c0 = pair1.second, c1 = pair2.second;
    double t_E = 0, t_L = 1;
    double dx = v1.x - v0.x, dy = v1.y - v0.y, dz = v1.z - v0.z;
    Color dc = Color(c1.r - c0.r, c1.g - c0.g, c1.b - c0.b);
    double x_min = -1, y_min = -1, z_min = -1;

    double x_max = 1, y_max = 1, z_max = 1;
    if (visible(dx, x_min-v0.x, t_E, t_L) && visible(-dx, v0.x-x_max, t_E, t_L)
        	&& visible(dy, y_min-v0.y, t_E, t_L) && visible(-dy, v0.y-y_max, t_E, t_L)
        		&& visible(dz, z_min-v0.z, t_E, t_L) && visible(-dz, v0.z-z_max, t_E, t_L)) {
        			isVisible = true;
					if (t_L < 1) {
						v1.x = v0.x + (dx * t_L);
						v1.y = v0.y + (dy * t_L);
						v1.z = v0.z + (dz * t_L);
						c1 = Color(dc.r * t_L + c0.r, dc.g * t_L + c0.g, dc.b * t_L+ c0.b);
					}
					if (t_E > 0) {
						v0.x = v0.x + (dx * t_E);
						v0.y = v0.y + (dy * t_E);
						v0.z = v0.z + (dz * t_E);
						c0 = Color(dc.r * t_E + c0.r, dc.g * t_E + c0.g, dc.b * t_E+ c0.b);
					}
				}

    pair1.first = v0;
    pair1.second = c0;
    pair2.first = v1;
    pair2.second = c1;

    return isVisible;
}

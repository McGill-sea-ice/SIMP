!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Module for the input grid domain use to                 !
!                    initialize the particule positions                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module domain

    !!!                      Variables definitions:                      !!!
    !   x_maxgeom : Upper limit of allowed x-regime                      [m]
    !   x_mingeom : Lower limit of allowed x-regime                      [m]
    !   y_maxgeom : Upper limit of allowed y-regime                      [m]
    !   y_mingeom : Lower limit of allowed y-regime                      [m]
    !   map_resolution : Side length of the map square grid              [m]
    !   x_numcell : number of cell in the x direction
    !   y_numcell : number of cell in the y direction
    implicit none

    !!! Note that the geometry variable are the coordinates at the       !!!
    !!! domain boundary cell center.                                     !!!

    double precision x_maxgeom, x_mingeom, y_maxgeom, y_mingeom ,          &
                     map_resolution
    parameter ( x_maxgeom = 435000.,                                           &
                x_mingeom = -195000.,                                          &
                y_maxgeom = 525000.,                                           &
                y_mingeom = -105000.,                                          &
                map_resolution = 30000)
    
    INTEGER :: x_numcell, y_numcell
    parameter ( x_numcell = 21,                                            &
                y_numcell = 21)

    INTEGER, dimension(x_numcell, y_numcell) :: map 
    
    parameter ( map = reshape((/   &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0, &
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  &
    /), (/x_numcell,y_numcell/)))
end module domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Module to be used for the nearest neighbour               !
!                       search based on tree structure                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tree_module

    use globals
    use data_variables

    implicit none
    
    !!!                        Local variables:                          !!!
    !   num_branch : Number of cell in each dimension
    !   area : Physical distance covered by the node
    !   leaves_id : ID of particle located on this branch
    !   child_node : Sublayer node linked to this branch 
    !   parent_node: Node from which the branch come
    !   child_branches : Array of branches linked to the node
    !   parent_branch : Branch from which the node come

    TYPE branch
        INTEGER, Dimension(:), Allocatable :: leaves_id(:)
        type(node), pointer :: child_node, parent_node
    end type branch

    type node
        INTEGER :: num_branch
        type(branch), pointer :: parent_branch, child_branches(:)
        double PRECISION :: area(dim,2)
    end type node


    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !              Subprograme to build the tree for the b-tree            !
    !                      nnps algorithm Cavelan 2020.                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine build_tree(n_particle, particle_id, domain, &
                                    current_node)

    use globals
    use data_variables

    implicit none
    
    !!!                        Tree structure :                          !!!
    !  Node layer:                              Root
    !  Branch layer:              /     /     |  ...  |     \     \
    !  Node layer:             node   node   node   node  node   node
    !  Branch layer:         / ... \    |   / ... \   |     |     |  \
    !  Node layer:        node     node   node
    !  Branch layer:     / | \     /...\    |
    !
    !    The leaves are located on the last branch in the tree structure    !

    !!!                        Local variables:                          !!!
    !   branch_size : Number of particle desired in each cell
    !   leaves_coor : Cell location of each particle
    !   alpha_t, beta_t : parameter to avoid over tree breadth explosion
    !   condition : Condition for while loop
    !   criteria : Criteria for rescaling n_branch
    !   i : loop index
    !   n_particle: number of particle in the node
    !   particle_id: ID of the particle in the node
    !   n_leaves : Number of particle in each cell/branche
    !   leaf_id : Array of particle ID that need to be add to branch type
    !   n_branch : Branching factor for the current node
    !   branch_id : ID of the child branches where a particle is located
    !   domain : current and next node domain
    !   min_hsml : Length to extend the domain so it encompasse particles
    !   current_node : Current layer where branches are generated
    !   next_node : Node associated with the new 
    !   new_branches: Array of pointer associated with the current node

    INTEGER          :: branch_size, n_particle, particle_id(n_particle),  &
                        leaves_coor(dim, n_particle), i, n_branch, branch_id
    INTEGER, ALLOCATABLE  :: n_leaves(:), leaf_id(:,:)
    DOUBLE PRECISION :: alpha_t, beta_t, criteria, domain(dim,2), min_hsml
    parameter         ( branch_size = 8,                                   &
                        alpha_t = 0.5,                                     &
                        beta_t  = 0.5 )
    LOGICAL          :: condition
    TYPE (node), pointer :: current_node, next_node
    TYPE (branch), pointer :: new_branches(:)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialize
    min_hsml = 0.
    leaves_coor = 0
    n_branch = 0
    criteria = 0.

    !Set while loop
    condition = .TRUE. 

    !Get the branching factor
    n_branch = ceiling((REAL(n_particle) / branch_size) ** (1. / dim) )

    !Update the node information:
    current_node%area = domain
    current_node%num_branch = n_branch

    !Create a layer of branches
    DO WHILE (condition)

        !Initialize criteria
        criteria = 0

        !Initialize n_leaves, leaf_id, etc.
        allocate( n_leaves(0 : n_branch ** dim - 1) )
        allocate( leaf_id(n_particle, 0 : n_branch ** dim - 1))
        n_leaves = 0
        leaf_id = 0
        criteria = 0

        !Loop on all particle in the current branch
        DO i = 1, n_particle

            !Get branch cooridnate for this particle
            leaves_coor(:,i) = floor( n_branch * (x(:,particle_id(i)) -    &
                               domain(:,1)) / (domain(:,2) - domain(:,1)) )
            
            !Add this particle to branch
            branch_id = leaves_coor(1,i) + leaves_coor(2,i) * n_branch
            n_leaves(branch_id) = n_leaves(branch_id) + 1

            !Store the particle ID in the branch
            leaf_id(n_leaves(branch_id), branch_id) = particle_id(i)

        ENDDO

        !Verify if the branching factor is small enough
        DO i = 0, n_branch ** 2 - 1
            IF (n_leaves(i) .LE. (alpha_t * branch_size) ) THEN
                criteria = criteria + 1 

            ENDIF
        ENDDO
        criteria = criteria / n_branch ** 2

        !Ajuste the branch factor if to big (can not be lower than 2)
        IF ( (criteria .GT. beta_t) .AND. (n_branch .GT. 2) ) THEN

            !Change branching factor
            n_branch = n_branch / 2
            
            !Make sure it is not lower than 2
            IF (n_branch .LT. 2) THEN
                n_branch = 2
            ENDIF

            !Update the current node information
            current_node%num_branch = n_branch

            !Deallocate and reinitialize varibales
            deallocate( n_leaves, leaf_id)
            leaves_coor = 0

        !If the branch factor is not to big then:
        ELSE

            !Allocate memory for the correct number of new branches
            allocate(new_branches(0 : n_branch ** dim - 1))

            !Create the branches and link them to the current node
            DO i = 0, n_branch ** 2 - 1

                !Initialize the current branch
                new_branches(i)%child_node => null()
                new_branches(i)%parent_node => current_node

                !Allocate memory for the leaves
                allocate(new_branches(i)%leaves_id(n_leaves(i)))

                !Assign the leaves to the branch
                new_branches(i)%leaves_id = leaf_id(:n_leaves(i),i)

            ENDDO

            !Link the current node to the new branch
            current_node%child_branches => new_branches

            !Update condition to stop
            condition = .FALSE.

        ENDIF
    ENDDO 

    !Clear undesired memory
    deallocate( n_leaves, leaf_id)

    !Create the new sub layer of node for each branch
    !Loop on branch
    DO i = 0, n_branch ** 2 - 1

        !Create a new node for cell with too much particle
        IF (size(new_branches(i)%leaves_id) .GT. branch_size) THEN

            !Create the new node 
            allocate(next_node)

            !Link the current branch and the nodes
            new_branches(i)%child_node => next_node
            next_node%parent_branch => new_branches(i)

            !Get domain of sublayer
            !The +/- hsml is used so domain encompasses all particles
            min_hsml = minval(hsml(particle_id))/100.
            domain = reshape(                                              &
                (/ MINVAL(x(1, new_branches(i)%leaves_id)) - min_hsml,     &
                   MINVAL(x(2, new_branches(i)%leaves_id)) - min_hsml,     &
                   MAXVAL(x(1, new_branches(i)%leaves_id)) + min_hsml,     &
                   MAXVAL(x(2, new_branches(i)%leaves_id)) + min_hsml /),  &
                   (/dim,2/) )

            !Redo the procedure for the subcell
            call build_tree(size(new_branches(i)%leaves_id),               &
                                 new_branches(i)%leaves_id, domain,        &
                                 next_node)

        ENDIF
    ENDDO

    end subroutine build_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Subprograme to set the current particle p and identify        !
    !            the possible interaction for this particle                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine p_interaction_tree(current_branch, switch, count_np)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particle to use                        [in]
    !   x      : Coordinates of all particles                           [in]
    !   hsml   : Smoothing lengths of particles                         [in]
    !   p_i  : List of first particle partner of interaction pair      [out]
    !   p_j  : List of second particle partner of interaction pair     [out]
    !   nip  : number of interacting particle                          [out]
    !   r_ij : distance between the interacting particles              [out]
    !   hsml_ij : mean smoothing length between particle i and j       [out]
    !   dx       : Distances in x and y between particles              [out]
    !   dvx      : Velocity difference in x and y between particles    [out]
    !   bdx     : Distances between particles and boundaries           [out]
    !   bp_i : List of boundary particle for boundary interaction      [out]
    !   bp_j : List of particle for boundary force interaction         [out]
    !   bnip : number of interacting particle for boundary force       [out]
    !   br_ij : distance between particles for boundary force          [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   ID_p : Current particle ID
    !   neighbor_range : Neighboring cell ID range
    !   branch_ID : Neighbor branch ID of the current particle
    !   count_np : number of interaction per particles                 [out]
    !   switch : Switch to avoid counting interaction twice
    !   j,k,l : Loop index
    !   current_branch : Current ending branch of the tree
    !   neighbor_branch : Branch with leaves that can interact with ID_p
    !   top_nod : Node that encompass the whole interaction domain of ID_p

    INTEGER             :: j, k, l, switch(ntotal), count_np(maxn),  &
                           ID_p, neighbor_range(dim, 2), branch_id
    TYPE(branch) :: current_branch, neighbor_branch
    TYPE(node)   :: top_node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialize
    neighbor_range = 0
    ID_p = 0
    branch_id = 0

    !Loop over all the particle in this cell
    DO j = 1, size(current_branch%leaves_id)

        !Get current particle ID
        ID_p = current_branch%leaves_id(j)

        !Assign the top_node
        top_node = current_branch%parent_node

        !Make sure the top_node encompass the range of interaction of ID_p
        DO WHILE ( ((x(1,ID_p) - hsml(ID_p)) .LT. top_node%area(1,1)) .OR. &
                   ((x(2,ID_p) - hsml(ID_p)) .LT. top_node%area(2,1)) .OR. &
                   ((x(1,ID_p) + hsml(ID_p)) .GT. top_node%area(1,2)) .OR. &
                   ((x(2,ID_p) + hsml(ID_p)) .GT. top_node%area(2,2)) )

            !Back track in the tree if top node area is not big enough
                   top_node = top_node%parent_branch%parent_node
        ENDDO

        !Get the range of possible interacting cell
        !Lower limit:
        neighbor_range(:, 1) = floor(                                      &
                        ( x(:,ID_p) - hsml(ID_p) - top_node%area(:,1) ) /  &
                        ( top_node%area(:,2) - top_node%area(:,1) ) *      &
                          top_node%num_branch )

        !Upper limit:
        neighbor_range(:, 2) = floor(                                      &
                        ( x(:,ID_p) + hsml(ID_p) - top_node%area(:,1) ) /  &
                        ( top_node%area(:,2) - top_node%area(:,1) ) *      & 
                          top_node%num_branch)

        !Compute interaction with all particle in the neighbor cell
        DO  k = neighbor_range(1, 1), neighbor_range(1, 2)
            DO l = neighbor_range(2, 1), neighbor_range(2, 2)

                !Get the branch ID
                branch_id = k + l * top_node%num_branch

                !Find all the ending branch and compute interaction
                neighbor_branch = top_node%child_branches(branch_id)

                call q_interaction_tree( neighbor_branch, ID_p, switch,    &
                 count_np)

            ENDDO
        ENDDO

        !Turn off particle to avoid counting the interaction twice
        switch(ID_p) = 0.
                
    ENDDO

    end subroutine p_interaction_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            Subprograme to loop on all neighbor branches from         !
    !       the particle ID_p and call the interaction when necessary      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine q_interaction_tree(c_branch, ID_p, switch,        &
                                            count_np)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!

    implicit none
    
    !!!                        Local variables:                          !!!
    !   c_branch : The current branch neighbor branch of particle ID_p
    !   count_np : number of interaction per particles                 [out]
    !   switch : Switch to avoid counting interaction twice
    !   ID_p : Current particle
    !   ID_q : Neighbor particle of ID_p
    !   i : loop index

    INTEGER :: i, switch(ntotal), count_np(maxn), ID_p, ID_q
    TYPE (branch) :: c_branch

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Verify if this branch is empty
    IF (size(c_branch%leaves_id) .GT. 0) THEN

        !If it is not the last branch loop on all child branch
        IF (associated(c_branch%child_node)) THEN

            !Compute interaction only if the branch is in the support domain
            IF ( .NOT. (                                                   &
                        ((x(1,ID_p) - hsml(ID_p)) .GT.                     &
                                    c_branch%child_node%area(1,2)) .OR.    &
                        ((x(2,ID_p) - hsml(ID_p)) .GT.                     &
                                    c_branch%child_node%area(2,2)) .OR.    &
                        ((x(1,ID_p) + hsml(ID_p)) .LT.                     &
                                    c_branch%child_node%area(1,1)) .OR.    &
                        ((x(2,ID_p) + hsml(ID_p)) .LT.                     &
                                    c_branch%child_node%area(2,1)))        &
                ) THEN

                !Loop on all branch of the child node
                DO i = 0, size(c_branch%child_node%child_branches(:)) - 1

                    !Loop on all child branches
                    call q_interaction_tree (                              &
                        c_branch%child_node%child_branches(i), ID_p,       &
                        switch, count_np)

                ENDDO
            ENDIF

        !If it is the last branch, compute interaction for its leaves
        ELSE

            !Loop over all the particle in the neighbor cell
            DO i = 1, size(c_branch%leaves_id)

                !Get neighbor particle ID
                ID_q = c_branch%leaves_id(i)

                !Don't count particle turned off:
                IF ((switch(ID_q) .EQ. 0) .OR. (ID_q .EQ. ID_p)) THEN
                    CYCLE
                ENDIF
                
                call interact(ID_p, ID_q, count_np)

            ENDDO
        ENDIF
    ENDIF

    end subroutine q_interaction_tree


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                 Subprograme to compute the interaction               !
    !                    between particles ID_p and ID_q                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine interact(ID_p, ID_q, count_np)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   p_i  : List of first particle partner of interaction pair      [out]
    !   p_j  : List of second particle partner of interaction pair     [out]
    !   nip  : number of interacting particle                          [out]
    !   r_ij : distance between the interacting particles              [out]
    !   hsml_ij : mean smoothing length between particle i and j       [out]
    !   dx       : Distances in x and y between particles              [out]
    !   dvx      : Velocity difference in x and y between particles    [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   d : distance of particle
    !   mhsml : mean smoothing length
    !   ID_p : Current particle
    !   ID_q : Neighbor particle of ID_p
    !   count_np : number of interaction per particles                 [out]

    INTEGER             :: ID_p, ID_q, count_np(maxn)
    DOUBLE PRECISION    :: d, mhsml

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Loop on all the branches
    !Compute distance and interaction length
    d     = SQRT(SUM((x(:,ID_p) - x(:,ID_q)) ** 2))
    mhsml = (hsml(ID_p)+ hsml(ID_q))/2.

    !Verify if particle are in their support domain:
    IF (d .LE. mhsml) THEN

        !Store interacting particle information
        nip         = nip + 1
        p_i(nip)    = ID_p
        p_j(nip)    = ID_q
        r_ij(nip)    = d
        hsml_ij(nip) = mhsml
        count_np(ID_p)  = count_np(ID_p) + 1
        count_np(ID_q)  = count_np(ID_q) + 1

    ENDIF

    end subroutine interact

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !               Subprograme to loop on all the tree structure          !
    !                 and call the interaction for each particle           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine tree_loop(c_node, switch, count_np)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!

    implicit none
    
    !!!                        Local variables:                          !!!
    !   c_node : The current node that will be forward track
    !   c_branch : The child branch to go trought
    !   next_node : Deeper node level in the tree
    !   count_np : number of interaction per particles                 [out]
    !   switch : Switch to avoid counting interaction twice
    !   i : loop index

    INTEGER :: i, switch(ntotal), count_np(maxn)
    TYPE (node):: c_node, next_node
    TYPE (branch) :: c_branch

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Loop on all the branches
    DO i = 0, c_node%num_branch ** 2 - 1

        !Verify if this branch is empty
        IF (size(c_node%child_branches(i)%leaves_id) .EQ. 0) THEN
            CYCLE
        ENDIF

        !Verify if this branch is the last (last if no child_node)
        IF (associated(c_node%child_branches(i)%child_node)) THEN

            !If not the last, go deeper in the tree
            next_node = c_node%child_branches(i)%child_node

            !Redo this subroutine with the child_node
            call tree_loop(next_node, switch, count_np)

        !If it is the last compute interaction for its leaves
        ELSE

            !Associate this branch to the current branch
            c_branch = c_node%child_branches(i)

            !Call interaction
            call p_interaction_tree(c_branch, switch, count_np)

        ENDIF
    ENDDO

    end subroutine tree_loop

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !              Subprograme to clear the tree structure memory          !
    !                                 allocated                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine delete_tree(c_node)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!

    implicit none
    
    !!!                        Local variables:                          !!!
    !   c_node : The current node that will be forward track
    !   next_node : Deeper node level in the tree
    !   i : loop index

    INTEGER :: i
    TYPE (node), pointer :: c_node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Loop on all the branches
    DO i = 0, c_node%num_branch ** 2 - 1

        !Verify if this branch is the last (last if no child_node)
        IF (associated(c_node%child_branches(i)%child_node)) THEN

            !If not the last, go deeper in the tree and redo
            call delete_tree(c_node%child_branches(i)%child_node)

            !Delete the leaves on the current branch
            deallocate(c_node%child_branches(i)%leaves_id)

        !If it is the last node layer delete all branches
        ELSE

            !Erase leaves of this branch
            deallocate(c_node%child_branches(i)%leaves_id)
            
        ENDIF
    ENDDO

    !Once all leaves deleted, delete the branches and node
    deallocate(c_node%child_branches)
    deallocate(c_node)

    end subroutine delete_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Subprograme to set the current boundary particle ID_p and     !
    !                call the interaction for this particle                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine boundary_tree_loop(root, count_np)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nbtotal : Total number of boundary particle                     [in]
    !   xb      : Coordinates of boundary particle                      [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   ID_p : Current particle ID
    !   count_np : number of interaction per particles                 [out]
    !   i : Loop index
    !   root : Top of the tree

    INTEGER :: i, ID_p, count_np(maxn)
    TYPE(node)   :: root

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialize
    ID_p = 0

    !Loop over all the boundary particle 
    DO i = 1, ntotalb

        !Get current particle ID
        ID_p = i

        !Make sure there is possible interaction with this b_particle
        IF ( .NOT. (xb(1,ID_p) .GT. root%area(1,1)) .AND. &
                    (xb(2,ID_p) .GT. root%area(2,1)) .AND. &
                    (xb(1,ID_p) .LT. root%area(1,2)) .AND. &
                    (xb(2,ID_p) .LT. root%area(2,2)) ) THEN

            CYCLE
        
        ELSE

            call bp_interaction_tree(root, ID_p, count_np)

        ENDIF
    ENDDO
    
    end subroutine boundary_tree_loop


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            Subprograme to loop on all neighbor branches from         !
    !           the boundary particle ID_p and call the interaction        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine bp_interaction_tree(c_node, ID_p, count_np)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   xb     : Coordinates of boundary particles                      [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   count_np : number of interaction per particles                 [out]
    !   ID_p : Current particle
    !   ID_q : Neighbor particle of ID_p
    !   i, j : loop index
    !   max_c_hsml : Maximum current smoothing length 

    INTEGER :: i, j, count_np(maxn), ID_p, ID_q
    DOUBLE PRECISION :: max_c_hsml
    TYPE (node) :: c_node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Loop on all branches and select only the one with possible interaction
    DO i = 0, size(c_node%child_branches) - 1

        !Verify if the branch is empty
        IF (size(c_node%child_branches(i)%leaves_id) .GT. 0) THEN

            !If it is not the last branch loop on all child branch
            IF (associated(c_node%child_branches(i)%child_node)) THEN

                !Compute interaction if the branch is in the support domain
                max_c_hsml = MAXVAL(hsml(c_node%child_branches(i)%leaves_id))
                IF ( .NOT. (                                               &
                        ((xb(1,ID_p) - max_c_hsml) .GT.                    &
                        c_node%child_branches(i)%child_node%area(1,2)) .OR.&
                        ((xb(2,ID_p) - max_c_hsml) .GT.                    &
                        c_node%child_branches(i)%child_node%area(2,2)) .OR.&
                        ((xb(1,ID_p) + max_c_hsml) .LT.                    &
                        c_node%child_branches(i)%child_node%area(1,1)) .OR.&
                        ((xb(2,ID_p) + max_c_hsml) .LT.                    &
                        c_node%child_branches(i)%child_node%area(2,1)))    &
                    ) THEN

                    !Go in the sublayer
                    call bp_interaction_tree (                              &
                        c_node%child_branches(i)%child_node, ID_p, count_np)

                ENDIF

            !If it is the last branch, compute interaction for its leaves
            ELSE

                !Loop over all the particle in the neighbor cell
                DO j = 1, size(c_node%child_branches(i)%leaves_id)

                    !Get neighbor particle ID
                    ID_q = c_node%child_branches(i)%leaves_id(j)
                    
                    call b_interact(ID_p, ID_q, count_np)

                ENDDO
            ENDIF
        ENDIF
    ENDDO
    
    end subroutine bp_interaction_tree


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                 Subprograme to compute the interaction               !
    !          between boundary particle ID_p and ice particle ID_q        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_interact(ID_p, ID_q, count_np)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   bdx     : Distances between particles and boundaries           [out]
    !   bp_i : List of boundary particle for boundary interaction      [out]
    !   bp_j : List of particle for boundary force interaction         [out]
    !   bnip : number of interacting particle for boundary force       [out]
    !   br_ij : distance between particles for boundary force          [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   d : distance of particle
    !   ID_p : Current particle
    !   ID_q : Neighbor particle of ID_p
    !   count_np : number of interaction per particles                 [out]

    INTEGER             :: ID_p, ID_q, count_np(maxn)
    DOUBLE PRECISION    :: d

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Loop on all the branches
    !Compute distance and interaction length
    d = SQRT(SUM((xb(:,ID_p) - x(:,ID_q)) ** 2))

    !Verify if particle are in their support domain:
    IF (d .LE. b_length) THEN

        !Store interacting particle information
        bnip        = bnip + 1
        bp_i(bnip)  = ID_p
        bp_j(bnip)  = ID_q
        br_ij(bnip) = d
        count_np(ID_q)  = count_np(ID_q) + 1

    ENDIF

    end subroutine b_interact

end module tree_module
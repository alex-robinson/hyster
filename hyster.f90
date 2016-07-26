
module hyster 

    implicit none 

    integer,  parameter :: dp  = kind(1.0d0)


    type hyster_class
        integer :: nt 
        real(dp), allocatable :: time(:), var(:), forc(:)
        real(dp) :: dvar_dt, dforc_dt
    end type 

    private

contains



    function hyster_init(nt) result(hyst)

        integer, intent(IN) :: nt 
        type(hyster_class) :: hyst 
        
        ! (Re)initialize hyster vectors
        if (allocated(hyst%time)) deallocate(hyst%time)
        if (allocated(hyst%var))  deallocate(hyst%var)
        if (allocated(hyst%forc)) deallocate(hyst%forc)
        allocate(hyst%time(nt))
        allocate(hyst%var(nt))
        allocate(hyst%forc(nt))


        return 

    end function hyster_init 

end module hyster 

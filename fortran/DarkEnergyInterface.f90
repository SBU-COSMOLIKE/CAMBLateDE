module DarkEnergyInterface
    use precision
    use interpolation
    use classes
    implicit none

    private

    type, extends(TCambComponent) :: TDarkEnergyModel
        logical :: is_cosmological_constant = .false.
        integer :: num_perturb_equations = 0
        !VM BEGINS
        real(dl) :: w_lam = -1_dl !VM not be used in Casarini except to init the search for effective constant w
        real(dl) :: wa    = 0._dl !VM not be used in Casarini except to init the search for effective constant w
        real(dl) :: cs2_lam = 1_dl
        logical  :: no_perturbations = .false.
        !VM ENDS

        contains
        procedure :: Init
        procedure :: BackgroundDensityAndPressure
        procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
        procedure :: diff_rhopi_Add_Term
        procedure :: PerturbationInitial
        !VM BEGINS
        !procedure :: PerturbationEvolve
        !VM ENDS
        procedure :: PrintFeedback
        ! do not have to implement w_de or grho_de if BackgroundDensityAndPressure is inherited directly
        procedure :: w_de
        procedure :: grho_de
        procedure :: Effective_w_wa !Used as approximate values for non-linear corrections
    end type TDarkEnergyModel

    public TDarkEnergyModel
    
    contains

    function w_de(this, a)
        class(TDarkEnergyModel) :: this
        real(dl) :: w_de, al
        real(dl), intent(IN) :: a

        w_de = -1._dl
    end function w_de  ! equation of state of the PPF DE

    function grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
        class(TDarkEnergyModel) :: this
        real(dl) :: grho_de, al, fint
        real(dl), intent(IN) :: a

        grho_de =0._dl
    end function grho_de

    subroutine PrintFeedback(this, FeedbackLevel)
        class(TDarkEnergyModel) :: this
        integer, intent(in) :: FeedbackLevel
    end subroutine PrintFeedback

    subroutine Init(this, State)
        use classes
        class(TDarkEnergyModel), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
    end subroutine Init

    subroutine BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
        !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
        class(TDarkEnergyModel), intent(inout) :: this
        real(dl), intent(in) :: grhov, a
        real(dl), intent(out) :: grhov_t
        real(dl), optional, intent(out) :: w

        if (this%is_cosmological_constant) then
            grhov_t = grhov * a * a
            if (present(w)) w = -1_dl
        else
            ! Ensure a valid result
            if (a > 1e-10) then
                grhov_t = grhov * this%grho_de(a) / (a * a)
            else
                grhov_t = 0._dl
            end if
            if (present(w)) w = this%w_de(a)
        end if
    end subroutine BackgroundDensityAndPressure

    subroutine Effective_w_wa(this, w, wa)
        class(TDarkEnergyModel), intent(inout) :: this
        real(dl), intent(out) :: w, wa

        w = -1
        wa = 0
    end subroutine Effective_w_wa

    subroutine PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
        class(TDarkEnergyModel), intent(inout) :: this
        real(dl), intent(out) :: dgrhoe, dgqe
        real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
        real(dl), intent(in) :: ay(*)
        real(dl), intent(inout) :: ayprime(*)
        integer, intent(in) :: w_ix

        dgrhoe=0
        dgqe=0
    end subroutine PerturbedStressEnergy

    function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
        class(TDarkEnergyModel), intent(in) :: this
        real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &
            k, grhov_t, z, k2, yprime(:), y(:), Kf1
        integer, intent(in) :: w_ix
        real(dl) :: ppiedot

        ! Ensure, that the result is set, when the function is not implemented by
        ! subclasses
        ppiedot = 0._dl
    end function diff_rhopi_Add_Term

    subroutine PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
        class(TDarkEnergyModel), intent(in) :: this
        real(dl), intent(inout) :: ayprime(:)
        real(dl), intent(in) :: a,adotoa, k, z, y(:), w
        integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve

    subroutine PerturbationInitial(this, y, a, tau, k)
        class(TDarkEnergyModel), intent(in) :: this
        real(dl), intent(out) :: y(:)
        real(dl), intent(in) :: a, tau, k
        !Get intinitial values for perturbations at a (or tau)
        !For standard adiabatic perturbations can usually just set to zero to good accuracy

        y = 0
    end subroutine PerturbationInitial

end module DarkEnergyInterface

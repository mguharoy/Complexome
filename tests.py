from Complexome import regulation, Perturbation


def test_perturbation_add_inits():
    assert (
        Perturbation.UNINITIALIZED + Perturbation.UNINITIALIZED
        == Perturbation.UNINITIALIZED
    )


def test_perturbation_add_init_up():
    assert (
        Perturbation.UNINITIALIZED + Perturbation.UP_REGULATED
        == Perturbation.UP_REGULATED
    )


def test_perturbation_add_up_init():
    assert (
        Perturbation.UP_REGULATED + Perturbation.UNINITIALIZED
        == Perturbation.UP_REGULATED
    )


def test_perturbation_add_init_down():
    assert (
        Perturbation.UNINITIALIZED + Perturbation.DOWN_REGULATED
        == Perturbation.DOWN_REGULATED
    )


def test_perturbation_add_down_init():
    assert (
        Perturbation.DOWN_REGULATED + Perturbation.UNINITIALIZED
        == Perturbation.DOWN_REGULATED
    )


def test_perturbation_add_init_altered():
    assert (
        Perturbation.UNINITIALIZED + Perturbation.ALTERED_COMPLEX
        == Perturbation.ALTERED_COMPLEX
    )


def test_perturbation_add_altered_init():
    assert (
        Perturbation.ALTERED_COMPLEX + Perturbation.UNINITIALIZED
        == Perturbation.ALTERED_COMPLEX
    )


def test_perturbation_add_ups():
    assert (
        Perturbation.UP_REGULATED + Perturbation.UP_REGULATED
        == Perturbation.UP_REGULATED
    )


def test_perturbation_add_downs():
    assert (
        Perturbation.DOWN_REGULATED + Perturbation.DOWN_REGULATED
        == Perturbation.DOWN_REGULATED
    )


def test_perturbation_add_altereds():
    assert (
        Perturbation.ALTERED_COMPLEX + Perturbation.ALTERED_COMPLEX
        == Perturbation.ALTERED_COMPLEX
    )


def test_perturbation_add_up_down():
    assert (
        Perturbation.UP_REGULATED + Perturbation.DOWN_REGULATED
        == Perturbation.ALTERED_COMPLEX
    )


def test_perturbation_add_up_altered():
    assert (
        Perturbation.UP_REGULATED + Perturbation.ALTERED_COMPLEX
        == Perturbation.ALTERED_COMPLEX
    )


def test_perturbation_add_down_up():
    assert (
        Perturbation.DOWN_REGULATED + Perturbation.UP_REGULATED
        == Perturbation.ALTERED_COMPLEX
    )


def test_perturbation_add_down_altered():
    assert (
        Perturbation.DOWN_REGULATED + Perturbation.ALTERED_COMPLEX
        == Perturbation.ALTERED_COMPLEX
    )


def test_regulation_up_1():
    assert regulation(list(range(1, 5))) == Perturbation.UP_REGULATED

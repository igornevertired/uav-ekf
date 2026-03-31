import numpy as np


def test_identification_quality_commercial():
    # Регрессионный тест: в линейной зоне (|alpha| <= ALPHA_MAX_SAFE) EKF должен сходиться
    # к разумной точности за 60 секунд.
    from feb.simulation.simulator import IntegratedSimulator

    sim = IntegratedSimulator(
        dt=0.02,
        accuracy_class="commercial",
        use_navigation=True,
        use_visualization=False,
        use_nav_for_identification=True,
    )
    ok = sim.run(T=60.0, profile_type="enhanced")
    assert ok

    results = sim.get_results()
    ident = results["identification"]

    # В реальном времени точность будет хуже; здесь проверяем, что модель/фильтр не деградировали.
    assert ident["mean_error"] < 10.0
    assert ident["max_error"] < 15.0

    # Проверяем, что alpha действительно был в линейной зоне
    true_hist = np.array(results["history"]["true_state"])
    alpha_deg_max = np.rad2deg(np.abs(true_hist[:, 1])).max()
    assert alpha_deg_max <= 12.0 + 1e-6


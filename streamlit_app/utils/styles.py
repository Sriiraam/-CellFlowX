import streamlit as st


def apply_global_style():
    st.markdown(
        """
        <style>

        /* APP BACKGROUND */
        .stApp {
            background:
                radial-gradient(circle at top left, #f4efff 0%, transparent 30%),
                radial-gradient(circle at top right, #f8f4ff 0%, transparent 35%),
                linear-gradient(180deg, #fcfaff 0%, #f7f3fc 100%);
            color: #2d2338;
        }

        /* MAIN WIDTH */
        .block-container {
            max-width: 1450px;
            padding-top: 1.6rem;
            padding-bottom: 3rem;
        }

        /* SIDEBAR */
        section[data-testid="stSidebar"] {
            background: linear-gradient(
                180deg,
                #f4edff 0%,
                #faf7ff 100%
            );
            border-right: 1px solid #e6dcf4;
        }

        section[data-testid="stSidebar"] * {
            font-size: 15px;
        }

        /* TYPOGRAPHY */
        html, body, [class*="css"] {
            font-family:
                Inter,
                -apple-system,
                BlinkMacSystemFont,
                "Segoe UI",
                sans-serif;
        }

        h1 {
            font-size: 2.35rem !important;
            font-weight: 800 !important;
            color: #281c35 !important;
        }

        h2 {
            font-size: 1.55rem !important;
            font-weight: 750 !important;
            color: #352643 !important;
            margin-top: 1.3rem !important;
        }

        h3 {
            font-size: 1.15rem !important;
            color: #473657 !important;
        }

        p, li {
            font-size: 16px;
            line-height: 1.7;
        }

        /* HERO */
        .hero {
            position: relative;
            overflow: hidden;
            padding: 2.6rem 2.8rem;
            border-radius: 26px;
            background:
                linear-gradient(
                    135deg,
                    #3b245d 0%,
                    #654693 55%,
                    #8363aa 100%
                );
            color: white;
            margin-bottom: 1.8rem;
            box-shadow: 0 14px 38px rgba(77, 52, 105, 0.18);
        }

        .hero:after {
            content: "";
            position: absolute;
            width: 330px;
            height: 330px;
            border-radius: 50%;
            background: rgba(255,255,255,0.08);
            right: -90px;
            top: -140px;
        }

        .hero-eyebrow {
            font-size: 0.82rem;
            font-weight: 700;
            letter-spacing: 0.16em;
            text-transform: uppercase;
            opacity: 0.8;
            margin-bottom: 0.55rem;
        }

        .hero-title {
            font-size: 3rem;
            font-weight: 850;
            letter-spacing: -0.04em;
            margin-bottom: 0.4rem;
        }

        .hero-subtitle {
            font-size: 1.08rem;
            line-height: 1.7;
            max-width: 850px;
            opacity: 0.91;
        }

        /* METRIC CARDS */
        .metric-card {
            background: rgba(255,255,255,0.84);
            border: 1px solid #e7def1;
            border-radius: 20px;
            padding: 1.25rem 1.35rem;
            min-height: 130px;
            box-shadow: 0 8px 26px rgba(70, 46, 95, 0.06);
            transition: all 0.25s ease;
        }

        .metric-card:hover {
            transform: translateY(-3px);
            box-shadow: 0 14px 34px rgba(70, 46, 95, 0.11);
            border-color: #ccbce1;
        }

        .metric-label {
            font-size: 0.83rem;
            font-weight: 650;
            color: #776886;
            margin-bottom: 0.4rem;
        }

        .metric-value {
            font-size: 2rem;
            font-weight: 800;
            color: #4f3570;
            letter-spacing: -0.03em;
        }

        .metric-subtitle {
            font-size: 0.82rem;
            color: #91839f;
            margin-top: 0.3rem;
        }

        /* CONTENT CARDS */
        .glass-card {
            background: rgba(255,255,255,0.82);
            border: 1px solid #e7def1;
            border-radius: 22px;
            padding: 1.55rem 1.65rem;
            box-shadow: 0 8px 28px rgba(76, 51, 105, 0.055);
            height: 100%;
        }

        .card-title {
            font-size: 1.05rem;
            font-weight: 750;
            color: #48315f;
            margin-bottom: 0.6rem;
        }

        .card-text {
            color: #665873;
            font-size: 0.95rem;
            line-height: 1.7;
        }

        /* SAMPLE CARDS */
        .sample-card {
            background: linear-gradient(
                145deg,
                rgba(255,255,255,0.94),
                rgba(247,241,255,0.94)
            );
            border: 1px solid #e4d8f2;
            border-radius: 20px;
            padding: 1.35rem;
            min-height: 175px;
            box-shadow: 0 8px 22px rgba(77, 52, 103, 0.06);
            transition: 0.25s ease;
        }

        .sample-card:hover {
            transform: translateY(-3px);
            box-shadow: 0 13px 30px rgba(77, 52, 103, 0.11);
        }

        .sample-id {
            font-size: 0.82rem;
            color: #8b7a99;
            font-weight: 650;
        }

        .sample-state {
            font-size: 1.12rem;
            font-weight: 800;
            color: #4f3570;
            margin-top: 0.45rem;
        }

        .sample-pct {
            font-size: 1.75rem;
            font-weight: 850;
            color: #72529a;
            margin-top: 0.7rem;
        }

        /* BUTTONS */
        .stButton > button {
            border-radius: 12px;
            border: 1px solid #d8c8e9;
            background: #f7f1ff;
            color: #4f3570;
            font-weight: 650;
            padding: 0.55rem 1rem;
        }

        .stButton > button:hover {
            border-color: #9f80bf;
            color: #3c2855;
        }

        /* EXPANDER */
        div[data-testid="stExpander"] {
            border: 1px solid #e4d9ef;
            border-radius: 15px;
            background: rgba(255,255,255,0.70);
        }

        /* INFO BOX */
        div[data-testid="stAlert"] {
            border-radius: 15px;
        }

        /* DIVIDER */
        hr {
            border-color: #e7dff0 !important;
        }

        </style>
        """,
        unsafe_allow_html=True,
    )


def metric_card(label, value, subtitle=""):
    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">{label}</div>
            <div class="metric-value">{value}</div>
            <div class="metric-subtitle">{subtitle}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

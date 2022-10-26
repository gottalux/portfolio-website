<template>
  <v-container fluid>
    <v-row class="text-left">
      <v-col cols="12" class="pa-0">
        <v-carousel
          cycle
          hide-delimiters
          next-icon="arrow_forward"
          prev-icon="arrow_back"
          :height="viewPortHeight"
        >
          <v-carousel-item
            v-for="(item, i) in art"
            :key="i"
            style="background-color: #ffeeee"
          >
            <v-img :src="item.picture" contain :height="viewPortHeight" />

            <div class="art-description">
              <p>{{ item.title }}</p>
            </div>
          </v-carousel-item>
        </v-carousel>
      </v-col>
    </v-row>
  </v-container>
</template>
<script>
import Vue from "vue";

export default {
  name: "ResearchPanel",
  data() {
    return {
      art: [
        {
          picture: require("../assets/aChaoticPresentation.png"),
          title: "A Chaotic Presentation",
        },
        {
          picture: require("../assets/backFromTheDead.png"),
          title: "Back From The Dead",
        },
        {
          picture: require("../assets/burnedOutGlow.jpg"),
          title: "Burned Out Glow",
        },
        {
          picture: require("../assets/dackel.png"),
          title: "Dackel",
        },
        {
          picture: require("../assets/dogAtTheEndOfTheRainbow.png"),
          title: "Dog at the End of the Rainbow",
        },
        {
          picture: require("../assets/dontEatMe.png"),
          title: "Don't Eat Me",
        },
        {
          picture: require("../assets/dumpling.png"),
          title: "Dumpling",
        },
        {
          picture: require("../assets/fingersAndOcToes.png"),
          title: "Fingers and Oc-Toes",
        },
        {
          picture: require("../assets/forceOfNature.png"),
          title: "Force of Nature",
        },
        {
          picture: require("../assets/frosty.png"),
          title: "Frosty",
        },
        {
          picture: require("../assets/mightyPirate.png"),
          title: "Mighty Pirate",
        },
        {
          picture: require("../assets/myUnbrokenHeart.png"),
          title: "My Unbroken Heart",
        },
        {
          picture: require("../assets/pigAtTheLavaPool.png"),
          title: "Pig at the Lava Pool",
        },
        {
          picture: require("../assets/spiltMilkDay.png"),
          title: "Spilt Milk (Day)",
        },
        {
          picture: require("../assets/spiltMilkNight.png"),
          title: "Spilt Milk (Night)",
        },
      ],
      windowHeight: window.innerHeight,
    };
  },
  computed: {
    viewPortHeight() {
      // The tabs at the top are 48 pixels in height
      return this.windowHeight - 48;
    },
  },
  created() {
    Vue.set(this, "art", this.randomiseArrayOrder(this.art));
    window.addEventListener("resize", () => {
      this.windowHeight = window.innerHeight;
    });
  },
  methods: {
    goToUrl(url) {
      window.open(url, "_blank");
    },
    randomiseArrayOrder(array) {
      let currentIndex = array.length,
        randomIndex;

      // While there remain elements to shuffle.
      while (currentIndex != 0) {
        // Pick a remaining element.
        randomIndex = Math.floor(Math.random() * currentIndex);
        currentIndex--;

        // And swap it with the current element.
        [array[currentIndex], array[randomIndex]] = [
          array[randomIndex],
          array[currentIndex],
        ];
      }

      return array;
    },
  },
};
</script>
<style lang="scss" scoped>
@import "../scss/main.scss";

.art-description {
  background-color: black;
  color: $adahlia-white;
  display: flex;
  height: 48px;
  justify-content: center;
  opacity: 0.75;
  position: absolute;
  transform: translateY(-48px);
  width: 100%;
}
</style>
